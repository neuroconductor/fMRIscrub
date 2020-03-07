#' Calculates PCA leverage or robust distance and identifies outliers.
#'
#' @param X A wide (observations x variables) numerical data matrix.
#' @param trend_filtering Should PCA Trend Filtering (PCATF) be used in place of
#'	regular PCA? Default is \code{TRUE}. If \code{TRUE}, \code{choosePCs} must
#'	be \{variance} and \code{method} must be \{leverage}.
#' @param trend_filtering.kwargs Options for \code{PCATF}, including the trend
#'	filtering parameter \code{lambda}.
#' @param choosePCs The method to use for choosing which PCs to retain:
#'	\code{variance} or \code{kurtosis}. Default is \code{kurtosis}.
#' @param kurt_quantile_cut The cutoff quantile to use if \code{choosePCs} is
#'	\code{kurtosis}. Default is .9.
#' @param kurt_detrend If \code{choosePCs} is \code{kurtosis}, should PCs be
#'	detrended before measuring kurtosis? Default is \code{TRUE}. Recommended if
#'	observations represent a time series.
#' @param method The method to use to measure outlyingness: \code{leverage},
#'	\code{robdist}, or \code{robdist_subset}. Default is \code{leverage}.
#' @param id_out If \code{TRUE} (default), will label outliers based on
#'	leverage or distance. Default is \code{TRUE}.
#' @param solve_directions Should the principal directions be solved for? These
#'	are needed to display the leverage images for outlying observations. Default
#'	is \code{TRUE}.
#' @param verbose Should occasional updates be printed? Default is \code{FALSE}.
#'
#' @return A clever object, i.e. a list with components
#' \describe{
#' 	\item{params}{A list with the \code{choosePCs}, \code{kurt_quantile_cut}, and \code{method} arguments used.}
#'  \item{PCs}{The PCs selected by \code{choosePCs}, i.e. the subsetted U matrix from the SVD.}
#'  \item{leverage}{The leverage of each observation. NULL if \code{method} is not PCA leverage.}
#'  \item{robdist}{The robust distance of each observation. NULL if \code{method} is not robust distance (subset).}
#'  \item{inMCD}{Whether each observation is within the MCD subset. NULL if \code{method} is not robust distance (subset).}
#'  \item{outliers}{An n X 3 data.frame indicating if each observation is an outlier at each of the three levels.}
#' }
#'
#' @importFrom robustbase rowMedians
#'
#' @export
#'
#' @examples
#' n_voxels = 1e4
#' n_timepoints = 100
#' X = matrix(rnorm(n_timepoints*n_voxels), ncol = n_voxels)
#'
#' clev = clever(X)
clever = function(
	X,
	trend_filtering = TRUE,
	trend_filtering.kwargs = NULL,
	choosePCs = c('variance', 'kurtosis'),
	kurt_quantile_cut = .9,
	kurt_detrend = TRUE,
	method = c('leverage','robdist_subset','robdist'),
	id_out = TRUE,
	solve_directions = TRUE,
	verbose = FALSE) {

	TOL <- 1e-8 # cutoff for zero variance/MAD detection

	# Return errors if args are not one of the acceptable options.
	choosePCs <- match.arg(choosePCs)
	method <- match.arg(method)

	if(choosePCs=='kurtosis'){
		if(trend_filtering){
			stop('A variance-based cutoff should be used with trend filtering.
				Set choosePCs=="variance" or trend_filtering=="FALSE"')
		}
		if(!is.numeric(kurt_quantile_cut)){
			stop('kurt_quantile_cut must be a number between 0 and 1.')
		}
		if((kurt_quantile_cut > 1) | (kurt_quantile_cut < 0)){
			stop('kurt_quantile_cut must be a number between 0 and 1.')
		}
	}

	choosePCs.fun <- switch(choosePCs, kurtosis=choosePCs_kurtosis,
																		 variance=choosePCs_variance,
																		 PCATF=choosePCs_PCATF)
	method.fun <- switch(method, leverage=PCleverage,
															 robdist=PCrobdist,
															 robdist_subset=PCrobdist_subset)
	id_out.fun <- switch(method, leverage=id_out.leverage,
															 robdist=id_out.robdist,
															 robdist_subset=id_out.robdist_subset)

	if(!is.matrix(X)){ X <- as.matrix(X) }
	N_ <- ncol(X)
	T_ <- nrow(X)
	if(N_ < T_) warning('Data matrix has more rows than columns.
		Check that observations are in rows and variables are in columns.')

	# Center and scale robustly.
	# Do it here directly instead of calling scale_med to save memory.
	if(verbose){ print('Centering and scaling.') }
	X <- t(X)
	## Center.
	X <- X - c(rowMedians(X, na.rm=TRUE))
	## Scale.
	mad <- 1.4826 * rowMedians(abs(X), na.rm=TRUE)
	zero_mad <- mad < TOL
	if(any(zero_mad)){
		if(all(zero_mad)){
			stop("All voxels are zero-variance.\n")
		} else {
			warning(cat("Warning: ", sum(zero_mad),
				" zero-variance voxels (out of ", length(zero_mad),
				"). These will be set to zero for estimation of the covariance.\n", sep=""))
			}
			mad[zero_mad] <- 1
	}
	X <- X/c(mad)
	X[zero_mad,] <- 0
	X <- t(X)
	rm(mad, zero_mad)

	# Compute the PC scores.
	if(verbose){ print('Computing the PC scores.') }
	if(trend_filtering){
		X.svd <- do.call(PCATF, c(list(X=X, K='mean var',
															solve_directions=solve_directions),
															trend_filtering.kwargs))
		rm(X)
		# Remove constant PCs.
		zero_var <- apply(X.svd$u, 2, var) < TOL
		if(any(zero_var)){
			print(paste0('PCATF returned ', sum(zero_var),
				' zero-variance PCs. Removing these.'))
			X.svd$u <- X.svd$u[,!zero_var]
			X.svd$d <- X.svd$d[!zero_var]
			X.svd$v <- X.svd$v[,!zero_var]
		}
	} else {
		if(solve_directions){
			X.svd <- svd(X)
			rm(X)
		} else {
			# Avoid computing U.
			XXt <- tcrossprod(X)
			rm(X)
			X.svd <- svd(XXt)
			X.svd$d <- sqrt(X.svd$d)
			X.svd$v <- NULL
			rm(XXt)
		}
	}
	gc()

	if(verbose){ print('Choosing PCs.') }
	# Choose which PCs to retain.
	#  First, get the keyword arguments...
	choosePCs.kwargs <- list(svd=X.svd)
	if((id_out == TRUE) & (method %in% c('robdist','robdist_subset'))){
		# Let q = N_/T_ (ncol(U)/nrow(U)). robustbase::covMcd()
		#  requires q <= approx. 1/2 for computation of the MCD covariance estimate.
		#  Higher q will use more components for estimation, thus retaining a
		#  higher resolution of information. Lower q will have higher breakdown
		#  points, thus being more resistant to outliers. (The maximal breakdown
		#  value is (N_ - T_ + 2)/2.) Here, we select q = 1/3 to yield a breakdown
		#  value of approx. 1/3. Since the subset method splits T_ into thirds, it
		#  must further reduce N_ by 1/3.
		q <- 1/3
		choosePCs.kwargs$max_keep <- max(1, floor(switch(method,
			robdist=T_*q,
			robdist_subset=T_*q/3)))
	}
	if(choosePCs == 'kurtosis'){
		choosePCs.kwargs$kurt_quantile_cut <- kurt_quantile_cut
		choosePCs.kwargs$detrend <- kurt_detrend
	}

	#  ...then call the respective functions for PC selection.
	if(trend_filtering){
		# We have already computed the PCs we want.
		if('max_keep' %in% names(choosePCs.kwargs)){
			chosen_PCs <- 1:min(ncol(X.svd$u), choosePCs.kwargs$max_keep)
		} else {
			chosen_PCs <- 1:ncol(X.svd$u)
		}
	} else {
		chosen_PCs <- do.call(choosePCs.fun, choosePCs.kwargs)
	}
	X.svd$u <- X.svd$u[,chosen_PCs]
	X.svd$d <- X.svd$d[chosen_PCs]
	if(!is.null(X.svd$v)){ X.svd$v <- X.svd$v[,chosen_PCs] }

	if(verbose){ print('Computing outlyingness.') }
	# Compute PCA leverage or robust distance.
	measure <- method.fun(X.svd$u)

	# Organize the output.
	PCs <- list(indices = chosen_PCs, svd=X.svd)
	params <- list(choosePCs=choosePCs, method=method,
								 trend_filtering=trend_filtering,
								 trend_filtering.kwargs=trend_filtering.kwargs)
	if(choosePCs == 'kurtosis'){
		params$kurt_quantile_cut = kurt_quantile_cut
		params$kurt_detrend = kurt_detrend
	}
	result <- list(params=params, PCs=PCs,
								 leverage=NULL, robdist=NULL, inMCD=NULL,
								 outliers=NULL, cutoffs=NULL)

	if(method == 'leverage'){
		result$leverage <- measure
	}
	if(method %in% c('robdist_subset','robdist')){
		inMCD <- measure$inMCD
		Fparam <- measure$Fparam
		measure <- measure$robdist
		result$robdist <- measure
		result$inMCD <- inMCD
	}

	if(id_out){
		if(verbose){ print('Identifying outliers.') }

		if(method == 'leverage'){
			id_out.kwargs <- list(leverage=measure)
		}
		if(method %in% c('robdist_subset','robdist')){
			id_out.kwargs <- list(distance=measure, inMCD=inMCD, Fparam=Fparam)
		}
		out <- do.call(id_out.fun, id_out.kwargs)
		result$outliers <- out$outliers
		result$cutoffs <- out$cutoffs
	}

	class(result) <- c('clever', class(result))

	return(result)
}
