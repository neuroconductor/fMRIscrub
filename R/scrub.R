#' Scrub fMRI data
#' 
#' Computes leverage, DVARS, or FD, and optionally flags the measure to identify
#'  artifactual time points.
#'
#' @param X If \code{method} is \code{"leverage"} or \code{"DVARS"}, this should
#'  be a \eqn{T} by \eqn{N} numeric matrix representing an fMRI run. There should
#'  not be any missing data (\code{NA} or \code{NaN}). Otherwise, if \code{method}
#'  is \code{"FD"}, this should be a \eqn{T} by \eqn{6} numeric matrix of the
#'  realignment parameters (3 translation and 3 rotation). 
#' @param method \code{"leverage"} (default), \code{"DVARS"} or \code{"FD"}
#' @param ... Additional arguments to each specific scrubbing function: 
#'  \code{\link{pscrub}}, \code{\link{DVARS}} or \code{\link{FD}}.
#' 
#' @return A list with components
#' \describe{
#'  \item{measure}{A length-T vector or data.frame with T rows, giving the outlyingness measure(s)}
#'  \item{measure_info}{Describes the outlyingness measure(s)}
#'  \item{outlier_cutoff}{The outlier cutoff value(s).}
#'  \item{outlier_flag}{A length-T vector or data.frame with T rows,  where \code{TRUE} indicates suspected outlier presence.}
#' }
#' 
#' @export
#' 
scrub <- function(X, method=c("leverage", "DVARS", "FD"), ...) {
  method <- match.arg(method, c("leverage", "DVARS", "FD"))
  FUN <- switch(method, leverage=pscrub, DVARS=DVARS, FD=FD)
  FUN(X, ...)
}

#' Scrub fMRI data in CIFTI format
#' 
#' Computes leverage, DVARS, or FD, and optionally flags the measure to identify
#'  artifactual time points. Requires the Connectime Workbench.
#' 
#' @param X Path to a CIFTI file, or a \code{"xifti"} object. 
#' @param method \code{"leverage"} or \code{"DVARS"}
#' @param brainstructures Character vector indicating which brain structure(s) 
#'  to use: \code{"left"} (left cortical surface), \code{"right"} (right 
#'  cortical surface) and/or \code{"subcortical"} (subcortical and cerebellar
#'  gray matter). Can also be \code{"all"} (obtain all three brain structures). 
#'  Default: \code{c("left", "right")} (excludes the subcortex). 
#' @param ... Additional arguments to each specific scrubbing function: 
#'  \code{\link{pscrub}} or \code{\link{DVARS}}.
#' 
#' @return A list with components
#' \describe{
#'  \item{measure}{A length-T vector or data.frame with T rows, giving the outlyingness measure(s)}
#'  \item{measure_info}{Describes the outlyingness measure(s)}
#'  \item{outlier_cutoff}{The outlier cutoff value(s).}
#'  \item{outlier_flag}{A length-T vector or data.frame with T rows,  where \code{TRUE} indicates suspected outlier presence.}
#' }
#' 
#' @export 
#' 
scrub_xifti <- function(X, method=c("leverage", "DVARS"), brainstructures=c("left", "right"), ...) {
  method <- match.arg(method, c("leverage", "DVARS"))
  FUN <- switch(method, leverage=pscrub, DVARS=DVARS)

  if (!requireNamespace("ciftiTools", quietly = TRUE)) {
    stop("Package \"ciftiTools\" needed. Please install it.", call. = FALSE)
  }

  if (!ciftiTools::is.xifti(X, messages=FALSE)) {
    is.fname <- length(X) == 1 && is.character(X) && file.exists(X)
    if (!is.fname) {stop("`X` must be a path to a CIFTI file or `xifti` object.")}
    X <- ciftiTools::read_xifti(X, brainstructures=brainstructures)
  } else {
    brainstructures <- match.arg(
      brainstructures,
      c("left","right","subcortical","all"),
      several.ok=TRUE
    )
    if ("all" %in% brainstructures) { 
      brainstructures <- c("left","right","subcortical")
    }
    if (!("left" %in% brainstructures)) { X <- ciftiTools::remove_xifti(X, "cortex_left") }
    if (!("right" %in% brainstructures)) { X <- ciftiTools::remove_xifti(X, "cortex_right") }
    if (!("subcortical" %in% brainstructures)) { X <- ciftiTools::remove_xifti(X, "subcortical") }
  }

  FUN(t(as.matrix(X)), ...)
}