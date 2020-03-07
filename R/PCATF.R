#' PCA Trend Filtering. From: https://github.com/Lei-D/PCATF
#'
#' @param X A numerical data matrix (observations x variables).
#' @param solve_directions Should the principal directions be solved for? These
#'	are needed to display the leverage images for outlying observations.
#' @param K The number of PCs to solve for, or \code{"mean var"}, in which case
#'	the number of PCs will be the amount of regular PCs whith variance above
#'	the mean.
#' @param lambda The trend filtering parameter; roughly, the filtering intensity.
#'	Default is .5 . Can be NULL (lets algorithm decide).
#' @param niter_max The number of iterations to use for approximating the PC.
#' @param tol The maximum 2-norm between iterations to accept as convergence.
#' @param verbose Print statements about convergence?
#'
#' @return SVD The trend-filtered SVD decomposition of X.
#'
#' @importFrom glmgen trendfilter
#' @importFrom far orthonormalization
#' @export
PCATF <- function(X, solve_directions = TRUE, K=NULL, lambda=.5,
									niter_max = 1000, tol = 1e-8, verbose=FALSE){

	N_ <- ncol(X)
	T_ <- nrow(X)
	X.svd_init <- svd(X)

	if(is.null(K)){
		K <- T_
	} else if(is.character(K)){
		if('mean var' %in% K){
			K <- choosePCs_variance(X.svd_init, max_keep=NULL, min_keep=NULL)
			K <- K[length(K)]
		}
	} else {
		if(!is.integer(K)){ print('K argument to PCATF is not recognized.') }
	}

	U <- matrix(NA, nrow = T_, ncol = K)
	D <- rep(NA, K)
	if(solve_directions){ V <- matrix(NA, nrow = N_, ncol = K) }

	for(k in 1:K){
		# Get initial eigenvector from regular svd.
		u <- X.svd_init$u[, k]
		d <- X.svd_init$d[k]
		v <- X.svd_init$v[, k]

		# Iterate between trendfiltering u and orthonormalizing v.
		for(i in 1:niter_max){
			u.last <- u
			tf.kwargs <- list(y = scale(X %*% v, center = FALSE, scale = d),
												x = 1:T_, nlambda = 1, k = 0)
			if(!is.null(lambda)){ tf.kwargs$lambda <- lambda }
			u <- do.call(glmgen::trendfilter, tf.kwargs)$beta
			v <- far::orthonormalization(
				crossprod(X, u), basis = FALSE, norm = TRUE)
			diff <- sqrt(mean((u - u.last)^2))
			if(diff < tol){
				break
			}
			if(verbose){
				if(i == niter_max){ print(paste0('PC ', k, ' did not converge.')) }
			}
		}

		d <- crossprod(u, X %*% v)[1, 1]
		X <- X - d * tcrossprod(u, v)
		U[, k] <- u
		D[k] <- d
		if(solve_directions){ V[, k] <- v}
	}
	out <- list(d = D, u = U)
	if(solve_directions){ out$v = V }
	return(out)
}
