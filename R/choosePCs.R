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
choosePCs_mean <- function(svd, max_keep=NULL, min_keep=NULL){
	U <- svd$u
	var <- svd$d

	# Identify how many PCs will be kept.
	n_keep <- sum(var > mean(var))

	# Constrain this number kept between the minumum and maximum, if specified.
	if(!is.null(max_keep)){n_keep <- min(n_keep, max_keep)}
	if(!is.null(min_keep)){n_keep <- max(n_keep, min_keep)}

	# PCs are already ordered by decreasing variance.
	U <- U[,1:n_keep]

	return(U)
}

#' Selects the principle components (PCs) of sufficient kurtosis from a SVD.
#'
#' First, the largest PCs which together explain 90% of the variance are
#' retained, and smaller components are removed. Each PC is detrended, and the
#' autocorrelation function of the top 10 PCs is measured. The kurtosis cutoff
#' is then the 90% quantile of the sampling distribution of kurtosis for
#' Normal data of the same length and autocorrelation as the PCs; it is
#' estimated by simulation or calculated from the theoretical asymptotic
#' distribution if the time series is long enough and autocorrelation is
#' negligible. Finally, the total number kept is constrained by the
#' \code{max_keep} and \code{min_keep} arguments.
#'
#' @param svd An SVD decomposition; i.e. a list containing u, d, and v.
#' @param max_keep If specified, the total number kept will be at most this
#' value.
#' @param min_keep If specified, the total number kept will be at least this
#' value.  Default 1.
#'
#' @return The subsetted u matrix with only the chosen columns (PCs).
#'
#' @importFrom e1071 kurtosis
#' @export
choosePCs_kurtosis <- function(svd, max_keep=NULL, min_keep=1, n_sim=5000){
	U <- svd$u
	m <- nrow(U)

	# First remove components that explain less than 90% of variation.
	cumvarexp <- cumsum(svd$d/sum(svd$d))
	n <- min(which((cumvarexp > .90)))
	U <- U[,1:n]

	# Detrend each PC.
	U.trend <- apply(U, 2, est_trend)
	U.detrended <- U - U.trend

	# Compute the kurtosis cutoff.
	U.top10 <- U.detrended[,1:min(10, n)]
	ACF <- est_ACF(U.top10, detrend=FALSE)
	ACF <- reg_ACF(ACF, method='AR')
	Sigma <- toeplitz(ACF)
	sim <- apply(sim_ts(n_sim, m, Sigma, fit_check=FALSE), 2, kurtosis, type=1)
	cut <- quantile(sim, .90)

	# Compute the kurtosis of remaining PCs.
	kurt <- apply(U.detrended, 2, kurtosis, type=1)

	# Identify how many PCs will be kept.
	n_keep <- sum(kurt > cut)

	# Constrain the number kept between the minumum and maximum, if specified.
	if(!is.null(max_keep)){n_keep <- min(n_keep, max_keep)}
	if(!is.null(min_keep)){n_keep <- max(n_keep, min_keep)}

	# The PCs with greatest kurtosis are chosen.
	to_keep <- order(-kurt)[1:n_keep]
	to_keep <- to_keep[order(to_keep)]
	U <- U[,to_keep]

	return(U)
}
