#' PCA Trend Filtering. From: https://github.com/Lei-D/PCATF
#'
#' @param X A numerical data matrix (observations x variables).
#' @param X.svd (Optional) The svd decomposition of X. Save time by providing
#'  this argument if the svd has already been computed. Default NULL.
#' @param solve_directions Should the principal directions be solved for? These
#'	will be needed to display the leverage images for outlying observations.
#' @param K (Optional) The number of PCs to solve for. If not provided, the
#'	number of PCs will be the amount of regular PCs with variance above
#'	the mean, up to 100 PCs. 
#' @param lambda The trend filtering parameter; roughly, the filtering intensity.
#'	Default is 0.5 . Can be NULL (lets algorithm decide).
#' @param niter_max The number of iterations to use for approximating the PC.
#' @param tol The maximum 2-norm between iterations to accept as convergence.
#' @param verbose Print statements about convergence?
#'
#' @return SVD The trend-filtered SVD decomposition of X (list with u, d, v).
#'
#' @importFrom glmgen trendfilter
#' @importFrom far orthonormalization
#' @export
PCATF <- function(X, X.svd=NULL, solve_directions = TRUE, K=NULL, lambda=.5,
  niter_max = 1000, tol = 1e-8, verbose=FALSE){

  stopifnot(is.numeric(X))
  if(is.null(X.svd)){ X.svd <- svd(X) }
  stopifnot(sort(names(X.svd))  == sort(c("u", "d", "v")))
  stopifnot(is.logical(solve_directions))
  if(is.null(K)){
    K <- length(choose_PCs.variance(X.svd, max_keep=NULL, min_keep=NULL))
    K <- min(100, K)
  } 
  stopifnot(is.numeric(K))
  stopifnot(K==round(K))
  stopifnot(is.numeric(lambda))
  if(lambda == 0){ 
    return(
      list(u = matrix(X.svd$u[, 1:K], ncol=K),
           d = X.svd$d[1:K],
           v = matrix(X.svd$v[, 1:K], ncol=K)
      )
    )
  }
  stopifnot(lambda > 0)
  stopifnot(is.numeric(niter_max))
  stopifnot(niter_max==round(niter_max))
  stopifnot(is.numeric(tol))
  stopifnot(tol > 0)
  stopifnot(is.logical(verbose))

  N_ <- ncol(X)
  T_ <- nrow(X)

  U <- matrix(NA, nrow = T_, ncol = K)
  D <- rep(NA, K)
  if(solve_directions){ V <- matrix(NA, nrow = N_, ncol = K) }

  for(k in 1:K){
    # Get initial eigenvector from regular svd.
    u <- X.svd$u[, k]
    d <- X.svd$d[k]
    v <- X.svd$v[, k]

    # Iterate between trendfiltering u and orthonormalizing v.
    for(i in 1:niter_max){
      u.last <- u
      tf.kwargs <- list(y = scale(X %*% v, center = FALSE, scale = d),
      x = 1:T_, nlambda = 1, k = 0)
      if(!is.null(lambda)){ tf.kwargs$lambda <- lambda }
      tryCatch({
        u <- do.call(glmgen::trendfilter, tf.kwargs)$beta
        v <- far::orthonormalization(crossprod(X, u), basis = FALSE, norm = TRUE)
      }, error = function(e) {
        # colinears vectors error
        u <- u.last
        if(verbose){print(paste0('Co-linear vectors: PC ', k, ' did not converge.'))
      }
      })
      diff <- sqrt(mean((u - u.last)^2))
      if(diff < tol){ break }
      if(verbose & i == niter_max){
        print(paste0('PC ', k, ' did not converge.'))
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
