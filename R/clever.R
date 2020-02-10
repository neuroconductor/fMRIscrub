#' Calculates PCA leverage or robust distance and identifies outliers.
#'
#' @param x A wide (obs x vars) data matrix of values.
#' @param choosePCs The method to use for choosing which PCs to retain. Default is kurtosis.
#' @param kurt_quantile_cut The cutoff quantile to use if \code{choosePCs} is \code{kurtosis}. Default is .9.
#' @param kurt_detrend If \code{choosePCs} is \code{kurtosis}, should PCs be detrended before measuring kurtosis?
#' 	Default is TRUE. Recommended if observations represent a time series.
#' @param method The method to use to measure outlyingness. Default is leverage.
#' @param id_out If TRUE (default), will label outliers based on leverage or distance.
#'
#' @return A clever object, i.e. a list with components
#' \describe{
#' 	\item{params}{A list with the \code{choosePCs}, \code{kurt_quantile_cut}, and \code{method} arguments used.}
#'  \item{PCs}{The PCs selected by \code{choosePCs}, i.e. the subsetted U matrix from the SVD.}
#'  \item{leverage}{The leverage of each observation. NULL if \code{method} is not PCA leverage.}
#'  \item{robdist}{The robust distance of each observation. NULL if \code{method} is not robust distance (subset).}
#'  \item{inMCD}{Whether each observation is within the MCD subset. NULL if \code{method} is not robust distance (subset).}
#'  \item{outliers}{An n x 3 data.frame indicating if each observation is an outlier at each of the three levels.}
#' }
#' @export
#'
#' @import stats
#' @importFrom robustbase covMcd
#'
#' @examples
#' n_voxels = 1e4
#' n_timepoints = 100
#' x = matrix(rnorm(n_timepoints*n_voxels), ncol = n_voxels)
#'
#' lev = clever(x)
clever = function(
	x,
	choosePCs = c('mean','kurtosis'),
	kurt_quantile_cut = .9,
	kurt_detrend = TRUE,
	method = c('leverage','robdist_subset','robdist'),
	id_out = TRUE) {

	choosePCs <- match.arg(choosePCs)  # return error if choosePCs arg not one of the acceptable options
	method <- match.arg(method)  # return error if method arg not one of the acceptable options

	if(choosePCs=='kurtosis'){
		if(!is.numeric(kurt_quantile_cut)){
			stop('kurt_quantile_cut must be a number between 0 and 1.')
		}
		if((kurt_quantile_cut > 1) | (kurt_quantile_cut < 0)){
			stop('kurt_quantile_cut must be a number between 0 and 1.')
		}
	}

	x <- as.matrix(x)
	p <- ncol(x)
	n <- nrow(x)
	if(p < n) warning('Data matrix has more rows than columns.
		Check that observations are in rows and variables are in columns.')

	# Center and scale robustly.
	x <- scale_med(x)

	# Perform dimension reduction.
	XXt <- (x %*% t(x))
	SVDi <- svd(XXt)

	# Choose which PCs to retain.
	choosePCs_kwargs <- list(svd=SVDi)
	choosePCs_fun <- switch(choosePCs, mean=choosePCs_mean, kurtosis=choosePCs_kurtosis)
	if((id_out == TRUE) & (method %in% c('robdist','robdist_subset'))){
		# Let q = n_PCs/n_timepoints (ncol(U)/nrow(U)). robustbase::covMcd()
		#  requires q <= approx. 1/2 for computation of the MCD covariance estimate.
		#  Higher q will use more components for estimation, thus retaining a
		#  higher resolution of information. Lower q will have higher breakdown
		#  points, thus being more resistant to outliers. (The maximal breakdown
		#  value is (n_PCs - n_timepoints + 2)/2.) Here, we select q = 1/3 to yield
		#  a breakdown value of approx. 1/3. Since the subset method splits
		#  n_timepoints into thirds, it must further reduce n_PCs by 1/3.
		q <- 1/3
		max_keep <- max(1, floor(switch(method,
			robdist=nrow(SVDi$u)*q,
			robdist_subset=nrow(SVDi$u)*q/3)))
		choosePCs_kwargs$max_keep <- max_keep
	}
	if(choosePCs == 'kurtosis'){
		choosePCs_kwargs$kurt_quantile_cut = kurt_quantile_cut
		choosePCs_kwargs$detrend = kurt_detrend
	}
	chosen_PCs <- do.call(choosePCs_fun, choosePCs_kwargs)
	Q <- ncol(chosen_PCs$U)

	# Compute PCA leverage or robust distance.
	method_fun <- switch(method, leverage=PCleverage,
		robdist_subset=PCrobdist_subset, robdist=PCrobdist)
	measure <- method_fun(chosen_PCs$U)

	# Organize the output.
	if(method %in% c('robdist_subset','robdist')){
		inMCD <- measure$inMCD
		Fparam <- measure$Fparam
		measure <- measure$robdist
	}
	params <- list(choosePCs=choosePCs, method=method)
	if(choosePCs == 'kurtosis'){ params$kurt_quantile_cut = kurt_quantile_cut }
	if(method == 'leverage'){
		result <- list(params=params, PCs=chosen_PCs,
			leverage=measure, robdist=NULL, inMCD=NULL, outliers=NULL, cutoffs=NULL)
	} else {
		result <- list(params=params, PCs=chosen_PCs,
			leverage=NULL, robdist=measure, inMCD=inMCD, outliers=NULL, cutoffs=NULL)
	}
	if(id_out){
		if(method=='leverage') id_out <- id_out.leverage(measure)
		if(method=='robdist_subset') id_out <- id_out.robdist_subset(measure, inMCD, Fparam)
		if(method=='robdist') id_out <- id_out.robdist(measure, inMCD, Fparam)
		result$outliers <- id_out$outliers
		result$cutoffs <- id_out$cutoffs
	}
	class(result) <- c('clever', class(result))

	return(result)
}
