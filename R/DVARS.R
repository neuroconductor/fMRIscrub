#' Estimates standard deviation robustly using the half IQR (and power trans.).
#'
#' @param x Numeric vector of data to estimate standard deviation for. 
#' @param d The scalar power transformation parameter. \eqn{w = x^{1/d}} is
#'  computed to obtain \eqn{w \sim N(\mu_w, \sigma_w^2)}.
#' @return Scalar for the robust estimate of standard deviation.
#' 
sd_hIQR <- function(x, d=1){
  w <- x^(1/d) # Power trans.: w~N(mu_w, sigma_w^2)
  sd <- (quantile(w, .5) - quantile(w, .25)) / (1.349/2) # hIQR
  out <- (d * median(w)^(d - 1)) * sd # Delta method
  # In the paper, the above formula incorrectly has d^2 instead of d.
  # The code on github correctly uses d.
  return(as.numeric(out))
}

#' Computes the DSE decomposition and DVARS-related statistics.
#' Citation: Insight and inference for DVARS (Afyouni and Nichols, 2018)
#' 
#' Differences from implementation at github.com/asoroosh/DVARS:
#' \itemize{
#'  \item The matrix is transposed.
#'  \item We center and scale the matrix differently (see \code{\link{scale_med}})
#'  \item We set all zero-variance voxels to zero during centering & scaling. This means that when we remove constant 0 or NA voxels, constant non-zero voxels are also removed.
#'  \item We use a tolerance of \eqn{1e-8} to detect non-zero voxels.
#' }
#'
#' @param X a T x N numeric matrix representing an fMRI run.
#' @param normalize Normalize the data as proposed in the original paper? Default is 
#'  \code{FALSE}.
#' @param norm_I The value to scale to. Default is \code{100}, as in the original
#'  paper.
#' @param verbose Should occasional updates be printed? Default is \code{FALSE}.
#'
#' @export
#' 
DVARS <- function(X, normalize=FALSE, norm_I=100, verbose=FALSE){
  T_ <- nrow(X)
  N_ <- ncol(X)

  if(normalize){
    # Normalization procedure from original DVARS paper and code.
    # Remove voxels of zeros (assume no NaNs or NAs)
    bad <- apply(X == 0, 2, all)
    if(any(bad)){
      if(verbose){ print(paste0('Zero voxels removed: ', sum(bad))) }
      X <- X[,!bad]
      N_ <- ncol(X)
    }
    
    # Scale the entire image so that the median average of each voxel is norm_I.
    X <- X / median(apply(X, 2, mean)) * norm_I

    # Center each voxel on its mean.
    X <- t(t(X) - apply(t(X), 1, mean))
  }

  # compute D/DVARS
  A_3D <- X^2
  Diff <- X[2:T_,] - X[1:(T_-1),]
  D_3D <- (Diff^2)/4
  A <- apply(A_3D, 1, mean)
  D <- apply(D_3D, 1, mean)
  DVARS_ <- 2*sqrt(D) # == sqrt(apply(Diff, 1, mean))

  # compute DPD
  DPD <- (D - median(D))/mean(A) * 100

  # compute z-scores based on X^2 dist.
  DV2 <- 4*D # == DVARS^2
  mu_0 <- median(DV2) # pg 305
  sigma_0 <- sd_hIQR(DV2, d=3) # pg 305: cube root power trans
  v <- 2*mu_0^2/sigma_0^2
  X <- v/mu_0 * DV2 # pg 298: ~X^2(v=2*mu_0^2/sigma_0^2)
  P <- pchisq(X, v)
  ZD <- ifelse(
    abs(P-.5)<.49999, # avoid overflow if P is near 0 or 1
    qnorm(1 - pchisq(X, v)), # I don't understand why they use 1-pchisq(X,v) instead of just pchisq(X,v)
    (DV2-mu_0)/sigma_0  # avoid overflow by approximating
  )

  out <- list(D=D, DVARS=DVARS_, DPD=DPD, ZD=ZD)
  return(out)
}
