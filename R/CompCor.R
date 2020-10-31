#' Hat matrix
#'
#' Get the hat matrix of X after adding a column of ones (for the intercept).
#'
#' @param X Numeric matrix. A column of ones will be added for the intercept.
#'
#' @return The hat matrix, a square matrix with the same number of rows/columns
#'  as the number of rows in \code{X}.
#'
#' @keywords internal
hat_mat <- function(X){
  X <- cbind(1, X)
  # The hat matrix is     X * (X^T * X)^(-1) * X^T.
  # This is identical to  X * solve(t(X) %*% X, t(X))
  # Could there be a more efficient way? 
  # https://stackoverflow.com/questions/9071020/compute-projection-hat-matrix-via-qr-factorization-svd-and-cholesky-factoriz/39298028#39298028
  X %*% solve(t(X) %*% X, t(X))
}

#' CompCor: get noise components
#'
#' @param X_noise The noise ROIs data
#' @param noise_nPC Number of PCs to obtain for each noise ROI
#'
#' @return A list with components X, X_noise, ROI_data, ROI_noise, noise_nPC,
#'  noise_erosion, noise_comps, and noise_var.
#' 
#' @keywords internal
CompCor.noise_comps <- function(X_noise, noise_nPC){

  N <- length(X_noise)
  noise_comps <- vector("list", N); names(noise_comps) <- names(X_noise)
  noise_var <- vector("list", N); names(noise_var) <- names(X_noise)

  for (ii in 1:N) {
    if (ncol(X_noise[[ii]]) == 0) { next }
    X_noise[[ii]] <- t(t(X_noise[[ii]]) - robustbase::rowMedians(t(X_noise[[ii]]), na.rm=TRUE))
    # Compute the PC scores.
    if (noise_nPC[[ii]] >= 1) {
      x <- svd(tcrossprod(X_noise[[ii]]), nu=noise_nPC[[ii]], nv=0)
      noise_comps[[ii]] <- x$u
      noise_var[[ii]] <- ((x$d^2)/sum(x$d^2))[1:noise_nPC[[ii]]]
    } else {
      x <- svd(tcrossprod(X_noise[[ii]]))
      noise_var[[ii]] <- (x$d^2)/sum(x$d^2)
      # Use enough PCs to explain the desired proportion of variance.
      noise_nPC[[ii]] <- min(length(x$d), sum(cumsum(noise_var[[ii]]) < noise_nPC[[ii]]) + 1)
      noise_comps[[ii]] <- x$u[,1:noise_nPC[[ii]],drop=FALSE]
      noise_var[[ii]] <- noise_var[[ii]][1:noise_nPC[[ii]]]
    }
  }

  list(noise_comps=noise_comps, noise_var=noise_var)
}

#' CompCor: regress
#'
#' @param X The data
#' @param noise_comps The noise components
#'
#' @return The data with the noise components regressed from it.
#' 
#' @keywords internal
CompCor.regress <- function(X, noise_comps){
  # Project each row of the data (which is why we transpose) on 
  #   the PCs (and the vector of ones for the intercept).
  noise_mat <- do.call(cbind, noise_comps)
  I_minus_H <- diag(nrow(noise_mat)) - hat_mat(noise_mat)
  t(I_minus_H %*% t(X))
}

#' CompCor
#'
#' The CompCor algorithm (Behzadi et. al., 2007) 10.1016/j.neuroimage.2007.04.042
#' 
#' The data are centered on each voxel timecourse's median.
#'
#' @inheritParams data_clever_CompCor_Params
#' @inheritParams noise_Params
#'
#' @return A list with entries \code{"data"} and \code{"noise"}
#'
#'  The entry \code{"data"} will be a \code{V x T} matrix where each row is a 
#'  data voxel (if it was originally an array, the voxels will be in spatial
#'  order) time series with each noise PC regressed from it. Otherwise, this 
#'  entry will be \code{NULL}.
#'
#'  The entry \code{"noise"} is a list of \code{T} \eqn{x} \code{noise_nPC} PC 
#'  scores, one for each \code{ROI_noise}.
#'
#' @importFrom robustbase rowMedians
#'
#' @export
CompCor <- function(X, ROI_data="infer", ROI_noise=NULL, noise_nPC=5, noise_erosion=NULL){

  out1 <- format_data(
    X=X, ROI_data=ROI_data, ROI_noise=ROI_noise, noise_nPC=noise_nPC, noise_erosion=noise_erosion
  )

  out2 <- CompCor.noise_comps(
    X_noise=out1$X_noise, noise_nPC=out1$noise_nPC
  )

  out1$X <- CompCor.regress(out1$X, out2$noise_comps)

  list(
    X=out1$X, 
    noise=list(PCs=out2$noise_comps, var=out2$noise_var, ROI=out1$ROI_noise)
  )
}