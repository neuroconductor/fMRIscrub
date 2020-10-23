#' Computes PCA leverage.
#'
#' Computes the leverage of each observation, in the PC score matrix (U).
#'  Optionally can identify the outliers.
#'
#' @param U An n x Q matrix of PC scores.
#' @param median_cutoff The outlier cutoff, in multiples of the median leverage.
#'  Default: \code{NA} (do not compute outliers).
#' 
#' @return A list with entries \code{"meas"} (the leverage values), 
#'  \code{"cut"} (the leverage cutoff value) and 
#'  \code{"flag"} (logical vector indicating the outliers). If 
#'  \code{!is.null(median_cutoff)}, all entries except \code{"meas"} are omitted.
#'  
#' @export
out_measures.leverage <- function(U, median_cutoff=NA){
  # Below line: same as diag(U %*% t(U)), but faster.
  lev <- apply(U, 1, function(x){sum(x*x)})
  out <- list(meas=lev)
  if (!is.null(median_cutoff)){
    out$cut <- median_cutoff * median(lev)
    out$flag <- out$meas > out$cut
  }
  out
}

#' Computes MCD distances.
#'
#' Computes robust minimum covariance determinant (MCD) distances across
#'  the observations (rows). The MCD method selects a subset of h observations
#'  whose covariance matrix has minimum determinant across all subsets of size
#'  h. The MCD distances are Mahalanobis distances using the estimates of
#'  center (mean) and scale (covariance matrix) based on that subset.
#'
#' @param U An n x Q matrix of PC scores.
#' @param quantile_cutoff The F-distribution quantile cutoff. Default: 
#'  \code{NA} (do not compute outliers).
#' 
#' @return A list with entries
#' \describe{
#'   \item{"meas"}{A vector of length n of with the robust distance estimate
#'    for each observation.}
#'   \item{"info"}{A list with entries "inMCD", "outMCD_scale", and "Fparam"}
#'   \item{"cut"}{The robust distance cutoff value}
#'   \item{"flag"}{Logical vector indicating the outliers}
#' }
#' 
#' If \code{is.null(quantile_cutoff)} the latter two elements are omitted.
#'
#' @importFrom MASS cov.mcd
#'
#' @export
out_measures.robdist <- function(U, quantile_cutoff=NA){ 
  n <- nrow(U)
  Q <- ncol(U)
  
  best <- c(cov.mcd(U)$best)
  inMCD <- 1:n %in% best
  U_in <- matrix(U[best,], ncol=Q) # observations that are used for the MCD estimates calculation
  xbar_star <- colMeans(U_in) # MCD estimate of mean
  U_ins <- scale(U_in, center = TRUE, scale = FALSE) 
  nU <- nrow(U_ins)
  S_star <- (t(U_ins) %*% U_ins)/(nU-1) # MCD estimate of covariance
  RD <- apply(U, 1, function(x) {t(x-xbar_star) %*% solve(S_star) %*% (x-xbar_star)})
  # Scale left-out observations to follow F-distribution.
  Fparam <- fit.F(Q, n, sum(inMCD))
  Fparam <- c(Fparam$c, Fparam$m, Fparam$df[1], Fparam$df[2])
  names(Fparam) <- c("c", "m", "df1", "df2")
  outMCD_scale <- Fparam["c"] * (Fparam["m"] - Q + 1) / (Q * Fparam["m"])   
  names(outMCD_scale) <- NULL
  
  out <- list(
    meas=RD, 
    info = list(inMCD=inMCD, outMCD_scale=outMCD_scale, Fparam=Fparam)
  )

  if (!is.null(quantile_cutoff)) {
    out$cut <- qf(p=quantile_cutoff, df1=Fparam["df1"], df2=Fparam["df2"])
    out$flag <- ifelse(inMCD, FALSE, RD * outMCD_scale > out$cut)
  }

  out
}