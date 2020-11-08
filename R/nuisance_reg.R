#' Nuisance regression
#'
#' Regresses matrix of nuisance variables from data matrix
#'
#' @param dat \eqn{N \times P} matrix of N measurements/timepoints, P data variables
#' @param nuis \eqn{N \times Q} matrix of N measurements/timepoints, Q nuisance regressors
#' @param intercept Add intercept column (constant) to \code{nuis}? Default: \code{TRUE}.
#' @return \code{dat} with \code{nuis} regressed from it
#' @export
nuisance_reg <- function(dat, nuis, intercept=TRUE){
  dat <- as.matrix(dat); nuis <- as.matrix(nuis)
  stopifnot(nrow(dat) == nrow(nuis))

  # Scale: I think this would help numeric stability.
  nuis <- scale(nuis)

  if (intercept) { nuis <- cbind(1/nrow(nuis), nuis) }

  I_mH <- diag(1200) - nuis %*% solve(t(nuis) %*% nuis, t(nuis))
  I_mH %*% dat
}