#' Estimates the trend of \code{ts} using a cubic spline.
#'
#' @param ts A numeric vector.
#' @param n_knots The number of knots to use for the spline model.
#'
#' @return The estimated trend.
#'
#' @import splines
#' @import stats
est_trend <- function(ts, n_knots=5){
	t = length(ts)
	i = 1:t

	loc_knots <- function(t, n_knots){
		if(n_knots >= t){stop('n_knots >= length(ts)')}
		return( t/(n_knots+1)*c(1:(n_knots)) )
	}

	# Do not use outliers (MAD > 3).
	out = abs(ts/stats::mad(ts)) > 3
	i[out] = NA

	# Use cubic spline with five knots.
	model = lm(ts~bs(i, knots=loc_knots(t, n_knots),
					 degree=3, Boundary.knots=c(1, t))
	)

	est = predict(model, newdata=data.frame(i=1:t))

	return(est)
}

#' Yields the upper diagonals of the input matrix.
#'
#' @param X The input matrix.
#' @param The diagonal index. i=1 is the main diagonal, i=2 is above it and
#'	shorter by 1, etc.
#'
#' @return The diagonal in a numerical vector.
upper_diag <- function(X,i) {
	n = nrow(X)
	len = n-i+1
	r = 1:len
	c = i:n
	idx = (c-1)*n+r
	return(X[idx])
}

#' Applies a summary function to each diagonal in X.
#'
#' @param X The input matrix.
#' @param FUN a function which takes in a numerical vector and returns a value.
#'
#' @return The summary function values in a vector, beginning with the one for
#' the main diagonal and ending with the one for the upper-right corner of X.
apply_diag <- function(X, FUN=mean) {
	sapply(1:nrow(X), function(i){FUN(upper_diag(X,i))})
}

# Cai, Ren, and Zhou, 2013
k_choice = function(beta, n, p){
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
	k = k_choice(beta, n, p)
	}
	k = min(k, n)
	# tapering
	if(type == 1){
		k = max(2, round(k/2)*2)
		w = c(rep(1, k/2), (k/2-1):0 / (k/2), rep(0, n-k))
	# banding
	} else if(type == 2){
		k = max(1, round(k))
		w = c(rep(1, k), rep(0, n-k))
	} else {
	stop('Invalid type.')
	}

	return(w)
}

#' Estimates the ACF of columns in X.
#'
#' @param X A txp matrix representing p length-t time series.
#' @param detrend Should the time series be detrended? Default TRUE.
#' @param cor_method The \code{method} parameter for \code{stats::cor} to use
#' 	for ACF estimation. Default is 'spearman'.
#' @param diag_center_method How ACF estimates for each column should be
#' 	aggregated. Default is \code{stats::median}.
#' @param beta The beta value for the tapering/banding estimator. See
#' 	\code{k_choice()}.
#' @param ACF_est_method Whether to use the 'tapering' or 'banding' estimator
#' 	from Cai, Ren and Zhou, 2013.
#' @param The k value for the tapering/banding estimator. Overrides \code{beta}.
#' 	See \code{k_choice()}.
#'
#' @return The ACF as a numeric vector.
est_ACF <- function(X,
	detrend=TRUE,
	cor_method='spearman', # pearson or spearman
	diag_center_method=median, #mean or median
	beta=.1,
	ACF_est_method='tapering',
	k=NULL){

	t = nrow(X)
	p = ncol(X)

	if(detrend){
		X_trend = apply(X, 2, est_trend)
		X = X - X_trend
	}

	# Obtain initial ACF estimate using a measure of center for each diagonal
	#  of the sample correlation matrix.
	acf = apply_diag(
		cor(t(X), method=cor_method),
		FUN=diag_center_method)

	# Stabilize higher-lag terms by tapering or banding.
	if(ACF_est_method == 'tapering'){
		w = acf_weights(p, t, beta, 1, k)
		acf = acf * w
	} else if(ACF_est_method == 'banding'){
		w = acf_weights(p, t, beta, 2, k)
		acf = acf * w
	} else if(ACF_est_method != 'none'){
		stop('ACF_est_method must be "tapering", "banding", or "none".')
	}

	return(acf)

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
sim_ts <- function(n, ts_length, Sigma=NULL,
                   fit_check=FALSE, fit_tol=.2, fit_tries=5, fit_print=FALSE){
  if(fit_tries < 1){ stop('Could not generate AR data within fit tolerance.') }

  if(is.null(Sigma)){Sigma = diag(ts_length)}

  # Generate X~N(0, I)
  X = t(mvrnorm(n, mu=rep(0, ts_length), diag(ts_length)))

  # Premultiply to obtain X~N(0, Sigma)
  if(!identical(Sigma, diag(ts_length))){
    V = Sigma #from AR_coefs: toeplitz(ARMAacf(ar=AR_coefs, lag.max=ts_length-1))
    V = svd(V)
    V_sqrt = V$u %*% diag(sqrt(V$d)) %*% t(V$v)
    X = V_sqrt %*% X
  }

  # Rudimentary check for proper AR model:
  # Try again if any (phi_est - phi_true) / phi_true > fit_tol. <-- fix: increase n
  if(fit_check){
    if(ts_length < 20){
      warning('Time series too short to estimate coefficients. Skipping check.')
    } else if(n < 20){
      warning('Too few runs to estimate coefficients. Skipping check.')
    } else {
      n_colsamp = min(1000, ncol(X))
      colsamp = sample(1:ncol(X), n_colsamp, replace=FALSE)
      AR_est = apply(X[,colsamp], 2, arima, order=c(length(AR_coefs), 0, 0),
                     include.mean=FALSE, method='CSS-ML')
      AR_est = apply(matrix(sapply(AR_est, '[[', 'coef'), ncol=n_colsamp), 1, median)
      if(fit_print){ print(list(phi_estimated=AR_est, phi_true=AR_coefs)) }
      errors = (AR_est - AR_coefs) / AR_coefs
      if(any(abs(errors) > fit_tol)){
        X = sim_ts(n=n, ts_length=ts_length, AR_coefs=AR_coefs,
                   fit_check=fit_check, fit_tol=fit_tol,
                   fit_tries=fit_tries-1, fit_print=fit_print)
      }
    }
  }
  return(X)
}
