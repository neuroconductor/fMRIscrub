#' Nuisance regression
#'
#' Regresses matrix of nuisance variables from data matrix
#'
#' @param dat \eqn{N \times P} matrix of N measurements/timepoints, P data variables
#' @param nuis \eqn{N \times Q} matrix of N measurements/timepoints, Q nuisance regressors.
#'  If \code{NULL}, only regress intercept (if \code{intercept==TRUE}. otherwise, return the original data).
#' @param intercept Add intercept column (constant) to \code{nuis}? Default: \code{TRUE}.
#' @return \code{dat} with \code{nuis} regressed from it
#' @export
nuisance_reg <- function(dat, nuis=NULL, intercept=TRUE){

  dat <- as.matrix(dat)

  if (is.null(nuis)) {
    if (!intercept) { return(dat) }
  } else {
    nuis <- as.matrix(nuis)
    stopifnot(nrow(dat) == nrow(nuis))
    # Scale: I think this would help numeric stability.
    nuis <- scale(nuis)
  }

  if (intercept) { nuis <- cbind(rep(1/nrow(dat), nrow(dat)), nuis) }

  I_mH <- diag(nrow(dat)) - nuis %*% solve(t(nuis) %*% nuis, t(nuis))
  I_mH %*% dat
}