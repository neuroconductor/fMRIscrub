#' Scrub fMRI data
#' 
#' Computes leverage, DVARS, or FD, and optionally flags the measure to identify
#'  artifactual time points.
#'
#' @param X a \eqn{T \times N} numeric matrix representing an fMRI run. There should
#'  not be any missing data (\code{NA} or \code{NaN}).
#' @param method \code{"leverage"} (default), \code{"DVARS"} or \code{"FD"}
#' @param ... Additional arguments to each specific scrubbing function: 
#'  \code{\link{clever}}, \code{\link{DVARS}} or \code{\link{FD}}.
#' 
#' @return A list with components
#' \describe{
#'  \item{measure}{A length-T vector or data.frame with T rows, giving the outlyingness measure(s)}
#'  \item{measure_info}{Describes the outlyingness measure(s)}
#'  \item{outlier_cutoff}{The outlier cutoff value(s).}
#'  \item{outlier_flag}{A length-T vector or data.frame with T rows,  where \code{TRUE} indicates suspected outlier presence.}
#' }
#' 
scrub <- function(X, method=c("leverage", "DVARS", "FD"), ...) {
  method <- match.arg(method, c("leverage", "DVARS", "FD"))
  FUN <- switch(method, leverage=clever, DVARS=DVARS, FD=FD)
  FUN(X, ...)
}