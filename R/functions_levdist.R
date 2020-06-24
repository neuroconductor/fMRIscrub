#' Computes PCA leverage.
#'
#' Computes the leverage of each observation, based on its PC scores.
#'
#' @param U An n x Q matrix of PC scores.
#'
#' @return A vector of length n with the leverage of each observation.
#' @export
PC.leverage <- function(U){
  return(diag(U %*% t(U)))
}

#' Computes MCD distances.
#'
#' Computes robust minimum covariance determinant (MCD) distances across
#'  the observations (rows).  The MCD method selects a subset of h observations
#'  whose covariance matrix has minimum determinant across all subsets of size
#'  h. The MCD distances are Mahalanobis distances using the estimates of
#'  center (mean) and scale (covariance matrix) based on that subset.
#'
#'
#' @param U An n x Q matrix of PC scores.
#'
#' @return A list with components
#' \describe{
#'   \item{robdist}{A vector of length n of with the robust distance estimate
#'    for each observation.}
#'  \item{inMCD}{A vector of length n indicating if each observation is within
#'    the MCD subset.}
#'  \item{Fparam}{The estimated parameters of the F distribution of MCD
#'    distances.}
#' }
#'
#' @importFrom MASS cov.mcd
#'
#' @export
PC.robdist <- function(U){
  n <- nrow(U)
  Q <- ncol(U)
  
  best <- c(cov.mcd(U)$best)
  inMCD <- 1:n %in% best
  U_in <- matrix(U[best,], ncol=Q) # observations that are used for the MCD estimates calculation
  xbar_star <- colMeans(U_in) # MCD estimate of mean
  U_ins <- scale(U_in, center = TRUE, scale = FALSE) # Damon: replace with U_ins - xbar_star?
  nU <- nrow(U_ins)
  S_star <- (t(U_ins) %*% U_ins)/(nU-1) # MCD estimate of covariance
  mah <- (apply(U, 1, function(x) t(x-xbar_star) %*% solve(S_star) %*% (x-xbar_star)))
  # Scale left-out observations to follow F-distribution.
  Fparam <- fit.F(Q, n, sum(inMCD))
  outMCD_scale <- Fparam$c * (Fparam$m - Q + 1) / (Q * Fparam$m) 
  #mah[!inMCD] <- outMCD_scale*mah[!inMCD] # Now, return scale and mah separately.
  
  result <- list(mah, inMCD, outMCD_scale, Fparam)
  names(result) <- c("mah", "inMCD", "outMCD_scale", "Fparam")
  return(result)
}