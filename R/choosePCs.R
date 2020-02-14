#' Selects the principle components of sufficient variance from a SVD.
#'

#' PCs with above-average variance are retained, but the total number kept is
#' constrained by the \code{max_keep} and \code{min_keep} arguments.
#'
#' @param svd An SVD decomposition; i.e. a list containing u, d, and v.
#' @param max_keep If specified, the total number kept will be at most this
#' value.
#' @param min_keep If specified, the total number kept will be at least this
#' value.
#'
#' @return The subsetted u matrix with only the chosen columns (PCs).
#' @export
choosePCs_variance <- function(svd, max_keep=NULL, min_keep=NULL){
	U <- svd$u
	var <- svd$d

	# Identify how many PCs will be kept.
	n_keep <- sum(var > mean(var))

	# Constrain this number kept between the minumum and maximum, if specified.
	if(!is.null(max_keep)){n_keep <- min(n_keep, max_keep)}
	if(!is.null(min_keep)){n_keep <- max(n_keep, min_keep)}

	# PCs are already ordered by decreasing variance.
	U <- U[,1:n_keep]

	return(list(U=U, indices=1:n_keep))
}

#' Selects the principle components (PCs) of sufficient kurtosis from a SVD.
#'
#' First, the largest PCs which together explain 90% of the variance are
#' retained, and smaller components are removed. Each PC is detrended. The
#' kurtosis cutoff is then the 90% quantile of the sampling distribution of
#' kurtosis for Normal data of the same length as the PCs; it is estimated by
#' simulation or calculated from the theoretical asymptotic distribution if the
#' time series is long enough. Finally, the total number kept is constrained by
#' the \code{max_keep} and \code{min_keep} arguments.
#'
#' @param svd An SVD decomposition; i.e. a list containing u, d, and v.
#' @param kurt_quantile_cut PCs with kurtosis of at least this quantile are kept.
#' @param detrend Should PCs be detrended before measuring kurtosis? Default is
#' 	TRUE. Recommended if observations represent a time series.
#' @param max_keep If specified, the total number kept will be at most this
#' value.
#' @param min_keep If specified, the total number kept will be at least this
#' value.  Default 1.
#' @param n_sim The number of simulation data to use for estimating the sampling
#' distribution of kurtosis.
#'
#' @return A list with the subsetted u matrix with only the chosen columns (PCs),
#' and the original indices of the PCs which were retained.
#'
#' @importFrom e1071 kurtosis
#' @importFrom MASS mvrnorm
#' @export
choosePCs_kurtosis <- function(svd, kurt_quantile_cut=.9, detrend=TRUE,
	max_keep=NULL, min_keep=1, n_sim=5000){
	U <- svd$u
	m <- nrow(U)

	# First remove components that explain less than 90% of variation.
	cumvarexp <- cumsum(svd$d/sum(svd$d))
	n <- min(which((cumvarexp > .90)))
	U <- U[,1:n]

	# Compute the kurtosis of the remaining PCs, detrending if applicable.
	if(detrend){
		U.dt <- U - apply(U, 2, est_trend)
		kurt <- apply(U.dt, 2, kurtosis, type=1)
	} else {
		kurt <- apply(U, 2, kurtosis, type=1)
	}

	# Determine the quantile cutoff.
	if(m < 1000){
		if(kurt_quantile_cut == .9){
			# Use precomputed empirical quantile.
			cut <- clever:::kurt_90_quant[m]
		} else {
			# Simulate and compute the empirical quantile.
			sim <- apply(t(mvrnorm(n_sim, mu=rep(0, m), diag(m))), 2, kurtosis, type=1)
			cut <- quantile(sim, kurt_quantile_cut)
		}
	} else {
		# Use theoretical quantile.
		cut <- qnorm(kurt_quantile_cut) * sqrt( (24*m*(m-1)^2) / ((m-3)*(m-2)*(m+3)*(m+5)) )
	}

	# Identify how many PCs will be kept.
	n_keep <- sum(kurt > cut)

	# Constrain the number kept between the minumum and maximum, if specified.
	if(!is.null(max_keep)){n_keep <- min(n_keep, max_keep)}
	if(!is.null(min_keep)){n_keep <- max(n_keep, min_keep)}

	# The PCs with greatest kurtosis are chosen.
	to_keep <- order(-kurt)[1:n_keep]
	to_keep <- to_keep[order(to_keep)]
	U <- U[,to_keep]

	return(list(U=U, indices=to_keep))
}
