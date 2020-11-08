#' Nuisance regression
#'
#' Regresses matrix of nuisance variables from data matrix
#'
#' @param data \eqn{N \times P} matrix of N measurements/timepoints, P data variables
#' @param nuis \eqn{N \times Q} matrix of N measurements/timepoints, Q nuisance regressors
#' @param intercept Add intercept column (constant) to \code{nuis}? Default: \code{TRUE}.
#' @return \code{data} with \code{nuis} regressed from it
#' @export
nuisance_reg <- function(data, nuis, intercept=TRUE){
  data <- as.matrix(data); nuis <- as.matrix(nuis)
  stopifnot(nrow(data) == nrow(nuis))
  
  if (intercept) { nuis <- cbind(1, nuis) }
  nuis <- apply(nuis, 2, scale) # I think this would help numeric stability.
  I_mH <- diag(1200) - nuis %*% solve(t(nuis) %*% nuis, t(nuis))
  I_mH %*% data
}