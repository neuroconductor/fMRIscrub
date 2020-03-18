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
#' @importFrom robustbase rowMedians
scale_med <- function(mat){
	TOL <- 1e-8
	# Transpose to use vector recycling (will revert after).
	mat <- t(mat)
	#	Center.
	mat <- mat - c(rowMedians(mat, na.rm=TRUE))
	# Scale.
	mad <- 1.4826 * rowMedians(abs(mat), na.rm=TRUE)
	zero_mad <- mad < TOL
	if(any(zero_mad)){
		if(all(zero_mad)){
			stop("All voxels are zero-variance.\n")
		} else {
			warning(paste0("Warning: ", sum(zero_mad),
				" zero-variance voxels (out of ", length(zero_mad),
				" ). These will be set to zero for estimation of the covariance.\n"))
		}
		mad[zero_mad] <- 1
	}
	mat <- mat/c(mad)
	mat[zero_mad,] <- 0
	# Revert transpose.
	mat <- t(mat)
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
#' @param ts A numeric vector to detrend.
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

#' Converts a vectorized matrix back to a volume time series.
Matrix_to_VolumeTimeSeries = function(mat, mask, sliced.dim = NA){
	in.mask = mask > 0
	t = nrow(mat)

	if(length(dim(mask)) == 3){
		dims = c(dim(mask), t)
	} else if(length(dim(mask)) == 2) {
		if(is.na(sliced.dim)){ sliced.dim=3 } #default to 3rd dim (axial)
		dims = switch(sliced.dim,
									c(1, dim(mask), t),
									c(dim(mask)[1], 1, dim(mask)[2], t),
									c(dim(mask), 1, t)
		)
	} else {
		stop('Not Implemented: mask must be 2D or 3D.')
	}

	vts = array(0, dim=dims)
	for(i in 1:t){
		vts[,,,i][in.mask] = mat[i,]
	}

	return(vts)
}
