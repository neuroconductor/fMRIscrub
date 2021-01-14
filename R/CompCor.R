#' CompCor: get noise components
#'
#' @param X_noise The noise ROIs data
#' @param center,scale,DCT,nuisance_too Center, scale, detrend, and nuisance
#'  regression
#' @param noise_nPC Number of PCs to obtain for each noise ROI
#'
#' @return A list with components X, X_noise, ROI_data, ROI_noise, noise_nPC,
#'  noise_erosion, noise_comps, and noise_var.
#' 
#' @keywords internal
CompCor.noise_comps <- function(X_noise, center, scale, DCT, nuisance_too, noise_nPC){
  TOL <- 1e-8

  N <- length(X_noise)
  noise_comps <- vector("list", N); names(noise_comps) <- names(X_noise)
  noise_var <- vector("list", N); names(noise_var) <- names(X_noise)
  noise_vartotal <- vector("list", N); names(noise_vartotal) <- names(X_noise)

  if (is.null(DCT)) { DCT <- 0 }
  detrend <- DCT > 0
  do_nreg <- !is.null(nuisance_too)
  if (do_nreg) {
    stopifnot(is.matrix(nuisance_too))
    stopifnot(nrow(nuisance_too) == T_)
  }

  if (length(noise_nPC) == 1) {
    noise_nPC <- as.list(rep(noise_nPC, N))
  } else {
    noise_nPC <- as.list(noise_nPC)
  }
  names(noise_nPC) <- names(X_noise)

  for (ii in 1:N) {
    T_ <- nrow(X_noise[[ii]])
    if (ncol(X_noise[[ii]])==0) { next }

    # Transpose.
    X_noise[[ii]] <- t(X_noise[[ii]])
    #	Center.
    if (center) { X_noise[[ii]] <- X_noise[[ii]] - c(rowMedians(X_noise[[ii]], na.rm=TRUE)) }
    # Detrend and perform nuisance regression.
    if (detrend | do_nreg) {
      B <- NULL
      if (detrend) { B <- dct_bases(T_, DCT) / sqrt((T_+1)/2) }
      if (do_nreg) { B <- cbind(B, nuisance_too) }
      if (center) {
        # Center design matrix robustly instead of using intercept term.
        B <- t(t(B) - c(rowMedians(t(B), na.rm=TRUE)))
      } else {
        B <- cbind(1, B)
      }
      X_noise[[ii]] <- t((diag(T_) - (B %*% t(B))) %*% t(X_noise[[ii]])) 
    }
    #	Center again for good measure.
    if (detrend && center) { X_noise[[ii]] <- X_noise[[ii]] - c(rowMedians(X_noise[[ii]], na.rm=TRUE)) }
    # Compute MADs.
    mad <- 1.4826 * rowMedians(abs(X_noise[[ii]]), na.rm=TRUE)
    X_constant <- mad < TOL
    if (any(X_constant)) {
      if (all(X_constant)) {
      stop(paste0("All data locations in noise ROI ", ii, " are zero-variance.\n"))
      } else {
        warning(paste0("Warning: ", sum(X_constant),
        " constant data locations (out of ", length(X_constant),
        ") in noise ROI ", ii, 
        ". These will be removed for estimation of the covariance.\n"))
      }
    }
    mad <- mad[!X_constant]; X_noise[[ii]] <- X_noise[[ii]][!X_constant,]
    # Scale.
    if (scale) { X_noise[[ii]] <- X_noise[[ii]]/c(mad) }
    # Revert transpose.
    X_noise[[ii]] <- t(X_noise[[ii]])
    if (ncol(X_noise[[ii]])==0) { next }

    # Compute the PC scores.
    x <- svd(tcrossprod(X_noise[[ii]]))
    noise_var[[ii]] <- x$d
    noise_vartotal[[ii]] <- sum(noise_var[[ii]])
    if (noise_nPC[[ii]] < 1) {
      # Use enough PCs to explain the desired proportion of variance.
      noise_nPC[[ii]] <- min(
        length(x$d), 
        sum(cumsum(noise_var[[ii]]) < noise_vartotal[[ii]]*noise_nPC[[ii]]) + 1
      )
    }
    noise_comps[[ii]] <- x$u[,seq(noise_nPC[[ii]]),drop=FALSE]
    noise_var[[ii]] <- noise_var[[ii]][seq(noise_nPC[[ii]])]
  }

  list(noise_comps=noise_comps, noise_var=noise_var, noise_vartotal=noise_vartotal)
}

#' CompCor: regress
#'
#' @param X The data
#' @param noise_comps The noise components
#'
#' @return The data with the noise components regressed from it.
#' 
#' @importFrom stats hat
#' @keywords internal
CompCor.regress <- function(X, noise_comps){
  # Project each row of the data (which is why we transpose) on 
  #   the PCs (and the vector of ones for the intercept).
  noise_mat <- do.call(cbind, noise_comps)
  I_minus_H <- diag(nrow(noise_mat)) - stats::hat(noise_mat)
  t(I_minus_H %*% t(X))
}

#' CompCor
#'
#' The CompCor algorithm (Behzadi et. al., 2007) 10.1016/j.neuroimage.2007.04.042
#' 
#' The data are centered on each voxel timecourse's median. If the data ROI is
#'  empty, the CompCor regressors will be calculated, but not regression will
#'  be performed.
#'
#' @inheritParams data_clever_CompCor_Params
#' @inheritParams noise_Params
#' @param center,scale Center the columns of the data by median, and scale the
#'  columns of the data by MAD? Default: \code{TRUE} for both. Affects both
#'  \code{X} and the noise data.
#' @param DCT Detrend the columns of the data using the discrete cosine
#'  transform (DCT)? Use an integer to indicate the number of cosine bases to 
#'  use for detrending. Use \code{0} (default) to forgo detrending. 
#' 
#'  The data must be centered, either before input or with \code{center}.
#' @param nuisance_too A matrix of nuisance signals to regress from the data
#'  before, i.e. a "design matrix." Should have \eqn{T} rows. Nuisance
#'  regression will be performed simultaneously with DCT detrending if 
#'  applicable. \code{NULL} to not add additional nuisance regressors. Affects 
#'  both \code{X} and the noise data.
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
CompCor <- function(
  X, ROI_data="infer", ROI_noise=NULL, 
  noise_nPC=5, noise_erosion=NULL,
  center=TRUE, scale=TRUE, DCT=0, nuisance_too=NULL
  ){

  out1 <- format_data(
    X=X, ROI_data=ROI_data, ROI_noise=ROI_noise, 
    noise_nPC=noise_nPC, noise_erosion=noise_erosion
  )

  out2 <- CompCor.noise_comps(
    X_noise=out1$X_noise, 
    center=center, scale=scale, DCT=DCT, nuisance_too=nuisance_too,
    noise_nPC=out1$noise_nPC
  )


  detrend <- DCT > 0
  do_nreg <- !is.null(nuisance_too)
  if (do_nreg) {
    stopifnot(is.matrix(nuisance_too))
    stopifnot(nrow(nuisance_too) == T_)
  }

  if (any(ROI_data)) {
    TOL <- 1e-8
    # Transpose.
    out1$X <- t(out1$X)
    #	Center.
    if (center) { out1$X <- out1$X - c(rowMedians(out1$X, na.rm=TRUE)) }
    # Detrend and perform nuisance regression.
    if (detrend | do_nreg) {
      B <- NULL
      if (detrend) { B <- dct_bases(T_, DCT) / sqrt((T_+1)/2) }
      if (do_nreg) { B <- cbind(B, nuisance_too) }
      if (center) {
        # Center design matrix robustly instead of using intercept term.
        B <- t(t(B) - c(rowMedians(t(B), na.rm=TRUE)))
      } else {
        B <- cbind(1, B)
      }
      out1$X <- t((diag(T_) - (B %*% t(B))) %*% t(out1$X)) 
    }
    #	Center again for good measure.
    if (center) { out1$X <- out1$X - c(rowMedians(out1$X, na.rm=TRUE)) }
    # Compute MADs.
    mad <- 1.4826 * rowMedians(abs(out1$X), na.rm=TRUE)
    X_constant <- mad < TOL
    if (any(X_constant)) {
      if (all(X_constant)) {
      stop(paste0("All data locations in noise ROI ", ii, " are zero-variance.\n"))
      } else {
        warning(paste0("Warning: ", sum(X_constant),
        " constant data locations (out of ", length(X_constant),
        ") in noise ROI ", ii, 
        ". These will be removed for estimation of the covariance.\n"))
      }
    }
    mad <- mad[!X_constant]; out1$X <- out1$X[!X_constant,]
    # Scale.
    if (scale) { out1$X <- out1$X/c(mad) }
    # Revert transpose.
    out1$X <- t(out1$X)

    out1$X <- CompCor.regress(out1$X, out2$noise_comps)
  } else {
    out1["X"] <- list(NULL)
  }

  list(
    X=out1$X, 
    noise=list(PCs=out2$noise_comps, var=out2$noise_var, ROI_noise=out1$ROI_noise)
  )
}