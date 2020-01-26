#' Estimates the trend of \code{ts} using a robust discrete cosine transform.
#'
#' @param ts A numeric vector.
#'
#' @param t The vector of values to regress.
#' @param robust Should a robust linear model be used? Default FALSE.
#'
#' @return The estimated trend.
#'
#' @importFrom robustbase lmrob
#' @importFrom robustbase lmrob.control
#' @export
est_trend <- function(ts, robust=TRUE){
  df <- data.frame(
    index=1:length(ts),
    ts=ts
  )

  i_scaled <- 2*(df$index-1)/(length(df$index)-1) - 1 #range on [-1, 1]

  df['p1'] <- cos(2*pi*(i_scaled/4 - .25)) #cosine on [-1/2, 0]*2*pi
  df['p2'] <- cos(2*pi*(i_scaled/2 - .5)) #cosine on [-1, 0]*2*pi
  df['p3'] <- cos(2*pi*(i_scaled*3/4  -.75)) # [-1.5, 0]*2*pi
  df['p4'] <- cos(2*pi*(i_scaled - 1)) # [2, 0]*2*pi

  if(robust){
		control <- lmrob.control(scale.tol=1e-3, refine.tol=1e-2) # increased tol.
		# later: warn.limit.reject=NULL
		trend <- lmrob(ts~p1+p2+p3+p4, df, control=control)$fitted.values
  } else {
		trend <- lm(ts~p1+p2+p3+p4, df)$fitted.values
	}

  return(trend)
}

#' Yields the upper diagonals of the input matrix.
#'
#' @param X The input matrix.
#' @param The diagonal index. i=1 is the main diagonal, i=2 is above it and
#'	shorter by 1, etc.
#'
#' @return The diagonal in a numerical vector.
upper_diag <- function(X,i) {
	n <- nrow(X)
	len <- n-i+1
	r <- 1:len
	c <- i:n
	idx <- (c-1)*n+r
	return(X[idx])
}

#' Applies a summary function to each upper diagonal in X.
#'
#' @param X The input matrix.
#' @param FUN a function which takes in a numerical vector and returns a value.
#'
#' @return The summary function values in a vector, beginning with the one for
#' the main diagonal and ending with the one for the upper-right corner of X.
apply_diag <- function(X, FUN=mean) {
	sapply(1:nrow(X), function(i){FUN(upper_diag(X,i))})
}

#' k from Cai, Ren, and Zhou, 2013
#'
#' @param beta The beta value (see paper).
#' @param n The number of observations of the data.
#' @param p The length of the data.
k_choice <- function(beta, n, p){
  return( (n*p/log(n*p))^(1/(2*beta+1)) )
}

#' Yields the tapering/banding weights from Cai, Ren and Zhou, 2013.
#'
#' If both k and beta are specified, k overrides beta.
#'
#' @param p The length of the data.
#' @param n The number of observations of the data.
#' @param beta The beta value (see paper).
#' @param type 1 for tapering, 2 for banding. Default 1.
#' @param k The k value (see paper). If specified, overrides \code{beta}.
#'
#' @return The weight vector.
acf_weights <- function(p, n, beta=.1, type=1, k=NULL){
	if(is.null(k)){
	k <- k_choice(beta, n, p)
	}
	k <- min(k, n)
	# tapering
	if(type == 1){
		k <- max(2, round(k/2)*2)
		w <- c(rep(1, k/2), (k/2-1):0 / (k/2), rep(0, n-k))
	# banding
	} else if(type == 2){
		k <- max(1, round(k))
		w <- c(rep(1, k), rep(0, n-k))
	} else {
	stop('Invalid type.')
	}

	return(w)
}

#' Regularizes an initial estimate of the ACF, \code{sigma}.
#'
#' @param method Regularization method: 'tapering' or 'banding' for the
#'	estimators from Cai, Ren and Zhou, 2013; or, 'AR' for an AR(6) model.
#' @param beta The beta value for the tapering/banding estimator. See
#' 	\code{k_choice()}.
#' @param k The k value for the tapering/banding estimator. Overrides \code{beta}.
#' 	See \code{k_choice()}.
#' @param p The number of time series used to estimate the input sigma. Needed for
#' 	tapering or banding estimator.
#'
#' @return The ACF as a numeric vector.
reg_ACF <- function(sigma, method, beta=.1, k=NULL, p=NULL){
	t <- length(sigma)

	# Stabilize higher-lag terms by tapering.
	if(method == 'tapering'){
		if(is.null(p)){stop('p must be specified for tapering estimator.')}
		w <- acf_weights(p, t, beta, 1, k)
		sigma <- sigma * w

	# Stabilize higher-lag terms by banding.
	} else if(method == 'banding'){
		if(is.null(p)){stop('p must be specified for banding estimator.')}
		w <- acf_weights(p, t, beta, 2, k)
		sigma <- sigma * w

	# Estimate an AR model, then obtain its ACF.
	} else if(method == 'AR'){
		ar <- solve(toeplitz(sigma[1:6]), sigma[2:7])
		sigma <- as.numeric(ARMAacf(ar=ar, lag.max=t-1))

	} else {
    stop(paste0('method ', as.character(method), ' is not recognized.'))
  }

	return(sigma)
}

#' Estimates the ACF of columns in X.
#'
#' @param X A numeric vector representing a length-t time series; or, a
#'	txp matrix representing p length-t time series.
#' @param detrend Should the time series be detrended? Default TRUE.
#' @param method 'matrixwise' will compute the multivariate correlation and
#'	take the center of the diagonals (default). 'vectorwise' will compute the
#'	ACF for each column, and then take the center at each lag. Different
#'	measures of correlation are available for 'matrixwise'; only Pearson
#'	correlation is available for 'vectorwise'.
#' @param cor_meas The \code{method} parameter for \code{stats::cor} to
#' 	use for ACF estimation with the 'matrixwise' method. Default is 'spearman'.
#' @param center_meas How values should be aggregated. Default is
#' \code{stats::median}.
#'
#' @return The ACF as a numeric vector.
est_ACF <- function(X,
	detrend=TRUE,
	method='matrixwise',
	cor_meas='spearman', # pearson or spearman, for multivariate X
	center_meas=median){ #mean or median, for multivariate X

	if(method=='vectorwise'){
		vector_ACF <- function(x){
			return(as.numeric(acf(x, lag.max=nrow(x)-1, plot=FALSE, demean=FALSE)$acf))
		}
		if(is.vector(X)){
			if(detrend){X <- X - est_trend(X)}
			sigma <- vector_ACF(X)
		} else {
			if(detrend){X <- X - apply(X, 2, est_trend)}
			sigma <- apply(X, 2, vector_ACF)
			sigma <- as.numeric(apply(sigma, 1, mean))
		}
	} else if(method=='matrixwise'){
		if(is.vector(X)){stop('Cannot compute matrixwise ACF on vector.')}
		if(detrend){ X <- X - apply(X, 2, est_trend) }
		# Obtain initial ACF estimate using a measure of center for each diagonal
		#  of the sample correlation matrix.
		sigma <- apply_diag(
			cor(t(X), method=cor_meas),
			FUN=center_meas
		)
	} else {
		stop(paste0(c('Unknown method argument in est_ACF(), ', method)))
	}
	return(sigma)
}

#' Simulates a time series with specified length and autocorrelation.
#'
#' @param n The number of time series to simulate.
#' @param ts_length The length of each time series.
#' @param Sigma The correlation matrix. Default is NULL, in which case the
#' 	identity matrix is used.
#' @param fit_check Should the time series be re-generated if the estimated
#' 	correlation deviates significantly from Sigma? Default FALSE.
#' @param fit_tol Try again if any (phi_est - phi_true) / phi_true > fit_tol,
#' 	with phi being the ACF coefficients. Only applies if fit_tol=TRUE.
#' @param fit_tries The maximum number of attempts. Only applies if fit_tol=TRUE.
#' @param fit_print Should the estimated ACF coefficients be displayed? (For
#' 	debugging.) Default FALSE. Only applies if fit_tol=TRUE.
#'
#' @return A ts_length by n matrix whose columns represent the simulated time
#' 	series.
#'
#' @importFrom MASS mvrnorm
sim_ts <- function(n, ts_length, Sigma=NULL,
                   fit_check=FALSE, fit_tol=.2, fit_tries=5, fit_print=FALSE){
  if(fit_tries < 1){ stop('Could not generate AR data within fit tolerance.') }

  if(is.null(Sigma)){Sigma <- diag(ts_length)}

  # Generate X~N(0, I)
  X <- t(mvrnorm(n, mu=rep(0, ts_length), diag(ts_length)))

  # Premultiply to obtain X~N(0, Sigma)
  if(!identical(Sigma, diag(ts_length))){
    V <- Sigma #from AR_coefs: toeplitz(ARMAacf(ar=AR_coefs, lag.max=ts_length-1))
    V <- svd(V)
    V_sqrt <- V$u %*% diag(sqrt(V$d)) %*% t(V$v)
    X <- V_sqrt %*% X
  }

  # Rudimentary check for proper AR model:
  # Try again if any (phi_est - phi_true) / phi_true > fit_tol. <-- fix: increase n
  if(fit_check){
	stop('Does not work at the moment.')
    if(ts_length < 20){
      warning('Time series too short to estimate coefficients. Skipping check.')
    } else if(n < 20){
      warning('Too few runs to estimate coefficients. Skipping check.')
    } else {
      n_colsamp <- min(1000, ncol(X))
      colsamp <- sample(1:ncol(X), n_colsamp, replace=FALSE)
      AR_est <- est_ACF(X, cor_meas='pearson', detrend=FALSE, center_meas=mean)
      if(fit_print){ print(list(phi_estimated=AR_est, phi_true=AR_coefs)) }
      errors <- (AR_est - AR_coefs) / AR_coefs
      if(any(abs(errors) > fit_tol)){
        X <- sim_ts(n=n, ts_length=ts_length, AR_coefs=AR_coefs,
                   fit_check=fit_check, fit_tol=fit_tol,
                   fit_tries=fit_tries-1, fit_print=fit_print)
      }
    }
  }
  return(X)
}
