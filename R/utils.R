#' Centers and scales a matrix robustly for the purpose of covariance estimation.
#'
#' Centers each column on its median, and scales each column by its median
#' absolute deviation (MAD). If any column MAD is zero, its values become zero
#' and a warning is raised. If all MADs are zero, an error is raised.
#'
#' @param mat A numerical matrix.
#'
#' @return The input matrix centered and scaled.
#'
#' @importFrom miscTools colMedians
scale_med <- function(mat){
	# mat is nxp; we want to scale the columns.
	n <- nrow(mat)
	p <- ncol(mat)

	# Center.
	mat <- sweep(mat, 2, colMedians(mat, na.rm=TRUE), '-')

	# Compute MAD and check for zero-variance voxels.
	mad <- 1.4826 * colMedians(abs(mat), na.rm=TRUE)
	zero_mad <- mad == 0
	if(any(zero_mad)){
		if(all(zero_mad)){
		stop("All voxels are zero-variance.\n")
		} else {
			warning(cat("Warning: ", sum(zero_mad),
				" zero-variance voxels (out of ", length(zero_mad),
				"). These will be set to zero for estimation of the covariance.\n", sep=""))
		}
	}

	# Scale.
	scale_col <- function(col, v){ return(ifelse(v != 0, col/v, 0)) }
	mat_scaled <- sweep(mat, 2, mad, scale_col)
	return(mat_scaled)
}

#' Computes the log likelihood of a sample of values from an F distribution.
#'
#' @param par A vector of length two which contains the degrees of freedom values.
#' @param vals A vector of values for which log likelihood will be calculated.
#' @param cutoff Values greater than the cutoff are removed before log likelihood calculation.
#'
#' @return A scalar which represents the log likelihood.
logL.F <- function(par, vals, cutoff){
	df1 <- par[1]
	df2 <- par[2]
	vals <- vals[vals <= cutoff]
	return(-1*sum(log(df(vals, df1, df2)) - log(pf(cutoff, df1, df2))))
}

#' Computes the log likelihood of a sample of values from a log normal distribution.
#'
#' @param par A vector of length two which contains the degrees of freedom values.
#' @param vals A vector of values for which log likelihood will be calculated.
#' @param cutoff Values greater than the cutoff are removed before log likelihood calculation.
#'
#' @return A scalar which represents the log likelihood.
logL.lnorm <- function(par, vals, cutoff){
	mean <- par[1]
	sd <- par[2]
	vals <- vals[vals <= cutoff]
	return(-1*sum(log(dlnorm(vals, mean, sd)) - log(plnorm(cutoff, mean, sd))))
}

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
