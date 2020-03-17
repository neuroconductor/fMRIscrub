#' Calculates PCA leverage or robust distance and identifies outliers.
#'
#' @param X A wide (observations x variables) numerical data matrix.
#' @param PCA_trend_filtering Should PCA Trend Filtering (PCATF) be used in place
#'	of regular PCA? Default is \code{TRUE}. If \code{TRUE}, \code{choose_PCs}
#'	must be \code{variance} and \code{method} must be \code{leverage}.
#' @param PCA_trend_filtering.kwargs Named list of arguments for \code{PCATF}:
#'	the trend filtering parameter \code{lambda} (Default 0.5), the number of
#'	iterations \code{niter_max} (Default 1000), convergence tolerance \code{tol}
#'	(Default 1e-8), and option to print updates \code{verbose} (Default FALSE).
#' @param choose_PCs The criteria for choosing which PCs to retain:
#'	\code{variance} or \code{kurtosis}. Default is \code{variance}. If trend
#'	filtering is being used, this must be \code{variance}.
#' @param kurt_quantile The cutoff quantile to use if \code{choose_PCs} is
#'	\code{kurtosis}. Default is 0.9.
#' @param kurt_detrend If \code{choose_PCs} is \code{kurtosis}, should PCs be
#'	detrended before measuring kurtosis? Default is \code{TRUE}. Recommended if
#'	observations represent a time series.
#' @param method The outlier measurement: \code{leverage}, \code{robdist},
#'	or \code{robdist_subset}. Default is \code{leverage}. If trend filtering
#'	is being used, this must be \code{leverage}.
#' @param id_out Should the outliers be identified? Default is \code{TRUE}.
#' @param solve_directions Should the principal directions be solved for? These
#'	are needed to display the leverage images for outlying observations. Default
#'	is \code{TRUE}. If the leverage images are not needed, this can be
#'	\code{FALSE} to reduce memory use.
#' @param verbose Should occasional updates be printed? Default is \code{FALSE}.
#'
#' @return A clever object, i.e. a list with components
#' \describe{
#' 	\item{params}{A list of all the arguments used.}
#'  \item{PCs}{
#'		\describe{
#'			\item{indices}{The indices of the selected PCs.}
#'			\item{svd}{The selected subset of the SVD.}
#'		}
#'	}
#'  \item{leverage}{The leverage of each observation. NULL if \code{method} is
#'		not PCA leverage.}
#'  \item{robdist}{The robust distance of each observation. NULL if
#'		\code{method} is not robust distance (subset).}
#'  \item{inMCD}{Whether each observation is within the MCD subset. NULL if
#'		\code{method} is not robust distance (subset).}
#'  \item{outliers}{An n X 3 data.frame indicating if each observation is an
#'		outlier at each of the three levels.}
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
	PCA_trend_filtering = TRUE,
	PCA_trend_filtering.kwargs = NULL,
	choose_PCs = c("variance", "kurtosis"),
	kurt_quantile = .9,
	kurt_detrend = TRUE,
	method = c("leverage","robdist_subset","robdist"),
	id_out = TRUE,
	solve_directions = TRUE,
	verbose = FALSE) {

	TOL <- 1e-8 # cutoff for detection of zero variance/MAD voxels

	# Argument Checks:
	choose_PCs <- match.arg(choose_PCs)
	method <- match.arg(method)
	if(!is.matrix(X)){ X <- as.matrix(X) }
	if(!is.logical(PCA_trend_filtering)){
		stop("Invalid argument: PCA_trend_filtering must be TRUE or FALSE.\n")
	}
	# 	If `PCA_trend_filtering`, check for variance cutoff and leverage method.
	if(PCA_trend_filtering){
		if(choose_PCs != "variance"){
			stop("Invalid argument: Trend Filtering requires choose_PCs==variance.\n")
		}
		if(method != "leverage"){
			stop("Invalid argument: Trend Filtering requires method==leverage.\n")
		}
	}
	#		If `choose_PCs==kurtosis`, check related arguments.
	if(choose_PCs=="kurtosis"){
		if(!is.numeric(kurt_quantile)){
			stop("Invalid argument: kurt_quantile must be numeric.\n")
		}
		if((kurt_quantile > 1) | (kurt_quantile < 0)){
			stop("Invalid argument: kurt_quantile must be between 0 and 1.\n")
		}
		if(!is.logical(kurt_detrend)){
			stop("Invalid argument: kurt_detrend must be TRUE or FALSE.\n")
		}
	}
	if(!is.logical(id_out)){
		stop("Invalid argument: id_out must be TRUE or FALSE.\n")
	}
	if(!is.logical(solve_directions)){
		stop("Invalid argument: solve_directions must be TRUE or FALSE.\n")
	}
	if(!is.logical(verbose)){
		stop("Invalid argument: verbose must be TRUE or FALSE.\n")
	}

	choose_PCs.fun <- switch(choose_PCs, kurtosis=choose_PCs.kurtosis,
																			 variance=choose_PCs.variance)
	method.fun <- switch(method, leverage=PC.leverage,
															 robdist=PC.robdist,
															 robdist_subset=PC.robdist_subset)
	id_out.fun <- switch(method, leverage=id_out.leverage,
															 robdist=id_out.robdist,
															 robdist_subset=id_out.robdist_subset)

	N_ <- ncol(X)
	T_ <- nrow(X)
	if(N_ < T_){
		warning("Data matrix has more rows than columns. Check that observations
			are in rows and variables are in columns.\n")
	}

	# Center and scale robustly.
	# Do it here instead of calling scale_med to save memory.
	if(verbose){ print("Centering and scaling the data matrix.") }
	X <- t(X)
	# Center.
	X <- X - c(rowMedians(X, na.rm=TRUE))
	# Scale.
	mad <- 1.4826 * rowMedians(abs(X), na.rm=TRUE)
	zero_mad <- mad < TOL
	if(any(zero_mad)){
		if(all(zero_mad)){
			stop("Error: All voxels are zero-variance. \n")
		} else {
			warning(paste0("Warning: ", sum(zero_mad),
				" zero-variance voxels (out of ", length(zero_mad),
				"). These will be set to zero for estimation of the covariance.\n"))
		}
		mad[zero_mad] <- 1
	}
	X <- X/c(mad)
	X[zero_mad,] <- 0
	X <- t(X)
	rm(mad, zero_mad)

	# Compute the PC scores.
	if(verbose){
		print(paste0("Computing the",
								 ifelse(PCA_trend_filtering, " trend-filtered", ""),
								 " PC scores",
							 	 ifelse(solve_directions, " and directions", ""), "."))
	}
	if(PCA_trend_filtering){
		X.svd <- do.call(PCATF, c(list(X=X, K="mean var",
																	 solve_directions=solve_directions),
															PCA_trend_filtering.kwargs))
		rm(X)
	} else {
		if(solve_directions){
			X.svd <- svd(X)
			rm(X)
		} else {
			# Avoid computing U.
			XXt <- tcrossprod(X)
			rm(X)
			X.svd <- svd(XXt)
			rm(XXt)
			X.svd$d <- sqrt(X.svd$d)
			X.svd$v <- NULL
		}
	}
	gc()

	# Remove constant PCs (likely present if trend filtering is used).
	zero_var <- apply(X.svd$u, 2, var) < TOL
	if(any(zero_var)){
		if(all(zero_var)){
			stop("Error: PCs are zero-variance.\n")
		}
		warning(paste0("Warning: ", sum(zero_var),
			" PCs are zero-variance. Removing these."))
		X.svd$u <- X.svd$u[,!zero_var]
		X.svd$d <- X.svd$d[!zero_var]
		X.svd$v <- X.svd$v[,!zero_var]
	}

	# Choose which PCs to retain.
	if(verbose){
		print(paste0("Choosing PCs with high ", choose_PCs, "."))
	}
	#  First, get the keyword arguments...
	choose_PCs.kwargs <- list(svd=X.svd)
	if((id_out == TRUE) & (method %in% c("robdist","robdist_subset"))){
		# Let q = N_/T_ (ncol(U)/nrow(U)). robustbase::covMcd()
		#  requires q <= approx. 1/2 for computation of the MCD covariance estimate.
		#  Higher q will use more components for estimation, thus retaining a
		#  higher resolution of information. Lower q will have higher breakdown
		#  points, thus being more resistant to outliers. (The maximal breakdown
		#  value is (N_ - T_ + 2)/2.) Here, we select q = 1/3 to yield a breakdown
		#  value of approx. 1/3. Since the subset method splits T_ into thirds, it
		#  must further reduce N_ by 1/3.
		q <- 1/3
		choose_PCs.kwargs$max_keep <- max(1, floor(switch(method,
			robdist=T_*q,
			robdist_subset=T_*q/3)))
	}
	if(choose_PCs == "kurtosis"){
		choose_PCs.kwargs$kurt_quantile <- kurt_quantile
		choose_PCs.kwargs$detrend <- kurt_detrend
	}
	#  ...then call the respective function for PC selection.
	if(PCA_trend_filtering){
		# We have already computed the PCs we want.
		if("max_keep" %in% names(choose_PCs.kwargs)){
			chosen_PCs <- 1:min(ncol(X.svd$u), choose_PCs.kwargs$max_keep)
		} else {
			chosen_PCs <- 1:ncol(X.svd$u)
		}
	} else {
		chosen_PCs <- do.call(choose_PCs.fun, choose_PCs.kwargs)
	}
	X.svd$u <- X.svd$u[,chosen_PCs]
	X.svd$d <- X.svd$d[chosen_PCs]
	if(!is.null(X.svd$v)){ X.svd$v <- X.svd$v[,chosen_PCs] }

	# Compute PCA leverage or robust distance.
	if(verbose){
		print(paste0("Computing ", method, " (the outlyingness measurement)."))
	}
	measure <- method.fun(X.svd$u)
	if(method %in% c("robdist_subset", "robdist")){
		inMCD <- measure$inMCD
		Fparam <- measure$Fparam
		measure <- measure$robdist
	}

	# Identify outliers.
	if(id_out){
		if(verbose){ print("Identifying outliers.") }
		if(method == "leverage"){
			id_out.kwargs <- list(leverage=measure)
		}
		if(method %in% c("robdist_subset", "robdist")){
			id_out.kwargs <- list(distance=measure, inMCD=inMCD, Fparam=Fparam)
		}
		out <- do.call(id_out.fun, id_out.kwargs)
	}

	# Organize the output.
	result <- list(params=NULL, PCs=NULL,
								 leverage=NULL, robdist=NULL, inMCD=NULL,
								 outliers=NULL, cutoffs=NULL)
	result$params <- list(PCA_trend_filtering=PCA_trend_filtering,
								 				PCA_trend_filtering.kwargs=PCA_trend_filtering.kwargs,
												choose_PCs=choose_PCs,
												kurt_quantile=kurt_quantile, kurt_detrend=kurt_detrend,
												method=method, id_out=id_out,
												solve_directions=solve_directions, verbose=verbose)
	result$PCs <- list(indices = chosen_PCs, svd=X.svd)
	if(method == "leverage"){
		result$leverage <- measure
	}
	if(method %in% c("robdist_subset", "robdist")){
		result$robdist <- measure
		result$inMCD <- inMCD
	}
	if(id_out){
		result$outliers <- out$outliers
		result$cutoffs <- out$cutoffs
	}
	class(result) <- c("clever", class(result))

	return(result)
}
