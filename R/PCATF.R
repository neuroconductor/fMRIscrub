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
#'	Default is 0.05 . Can be NULL (lets algorithm decide).
#' @param niter_max The number of iterations to use for approximating the PC.
#' @param TOL The maximum 2-norm between iterations to accept as convergence.
#' @param verbose Print statements about convergence?
#'
#' @return SVD The trend-filtered SVD decomposition of X (list with u, d, v).
#'
#' @importFrom glmgen trendfilter
#' @examples
#' set.seed(12345)
#' U = matrix(rnorm(100*3),ncol=3)
#' U[20:23,1] = U[20:23,1] + 3
#' U[40:43,2] = U[40:43,2] - 2
#' U = svd(U)$u
#' D = diag(c(10,5,1))
#' V = svd(matrix(rnorm(3*20),nrow=20))$u
#' X = U %*% D %*% t(V)
#' out3 = PCATF(X, K=3, lambda=.75)
#' matplot(out3$u, ty='l')
#' out3$d
#' plot(rowSums(out3$u^2), ty='l')
#'
#' # Orthonormalized
#' out3_svd = svd(out3$u)
#' matplot(out3_svd$u, ty='l')
#' out3_svd$d
#' plot(rowSums(out3_svd$u^2), ty='l')
#' @export
PCATF <- function(X, X.svd=NULL, solve_directions = TRUE, K=NULL, lambda=.05,
                  niter_max = 1000, TOL = 1e-8, verbose=FALSE){

  stopifnot(is.numeric(X))
  if(is.null(X.svd)){ X.svd <- svd(X) }
  stopifnot(sort(names(X.svd))  == sort(c("u", "d", "v")))
  stopifnot(is.logical(solve_directions))
  if(is.null(K)){
    K <- length(choose_PCs.variance(X.svd))
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
  stopifnot(is.numeric(TOL))
  stopifnot(TOL > 0)
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
      tf.kwargs <- list(y = X %*% v, #scale(X %*% v, center = FALSE, scale = d),
                        x = 1:T_, nlambda = 1, k = 0)
      if(!is.null(lambda)){ tf.kwargs$lambda <- lambda }
      u <- do.call(glmgen::trendfilter, tf.kwargs)$beta
      u <- u / sqrt(sum(u^2))
      if(any(is.na(u))){
        u <- u.last
        break
      }
      v <- crossprod(X, u)
      v <- v / sqrt(sum(v^2))
      diff <- sqrt(mean((u - u.last)^2))
      if(diff < TOL){ break }
      if(verbose & i == niter_max){
        cat(paste0('PC ', k, ' did not converge.\n'))
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
