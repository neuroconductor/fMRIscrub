#' Framewise Displacement
#'
#' Calculates Framewise Displacement (FD) 
#' 
#' The FD formula is taken from (Power, 2012):
#'
#'  \eqn{FD_i = \lvert \Delta x_i \rvert + \lvert \Delta y_i \rvert + \lvert \Delta z_i + \lvert \Delta \alpha_i \rvert + \lvert \Delta \beta_i \rvert + \lvert \Delta \gamma_i \rvert}, 
#'  where \eqn{i} is the timepoint; \eqn{x}, \eqn{y} and \eqn{z} are the 
#'  translational realignment parameters (RPs);
#'  \eqn{\alpha}, \eqn{\beta} and \eqn{\gamma} are the rotational RPs;
#'  and \eqn{\Delta x_i = x_{i-1} - x_i} (and similarly for the other RPs).
#'
#'  Citation: 10.1016/j.neuroimage.2011.10.018
#' 
#' @param X An \eqn{N \times 6} matrix in which the first three columns represent the
#'  translational RPs (\code{trans_units}), and the second three columns represent
#'  the rotational RPs (\code{rot_units}). If \code{rot_units} measures an angle,
#'  it will be converted to \code{trans_units} by measuring displacement on a 
#'  sphere of radius \code{brain_radius} \code{trans_units}.
#'
#'  Alternatively, this can be the file path to an \eqn{N \times 6} matrix which can be
#'  read with \code{\link[utils]{read.table}} (fields separated by white-space; no
#'  header).
#' @param trans_units \code{"mm"} for millimeters (default), \code{"cm"} 
#'  for centimeters, or \code{"in"} for inches.
#' @param rot_units \code{"deg"} for degrees (default), \code{"rad"} for radians,
#'  or one of the \code{trans_units} options.
#' @param brain_radius If \code{rot_units} measures an angle, the rotational RPs
#'  are transformed to a spatial measurement representing the displacement on a 
#'  sphere of radius \code{brain_radius} \code{trans_units}.
#' 
#'  If \code{brain_radius} is \code{NULL} (default), its value will be set to 
#'  (the equivalent of) 50 mm.
#' @param detrend Detrend each RP with \code{\link{dct_mat}} before computing FD?
#'  Default: \code{FALSE}.
#' @return A length \eqn{N} vector of FD values in \code{trans_units}.
#'
#' @importFrom utils read.table
#' @export
FD <- function(
  X, trans_units = c("mm", "cm", "in"), rot_units = c("deg", "rad", "mm", "cm", "in"), 
  brain_radius=NULL, detrend=FALSE) {

  if (is.character(X)) { X <- read.table(X) }
  X <- as.matrix(X); stopifnot(is.matrix(X))
  stopifnot(nrow(X) > 1); stopifnot(ncol(X) >= 6)
  if (ncol(X) > 6) { 
    warning(paste(
      "`X` has more than 6 columns.",
      "using the first 3 as translation RPs,",
      "the second 3 as rotation RPs, and discarding the rest.\n"
    ))
    X <- X[,1:6]
  }

  # Convert RPs and brain radius to mm.
  trans_units <- match.arg(trans_units, trans_units)
  X[,1:3] <- X[,1:3] * switch(trans_units, mm=1, cm=10, `in`=25.4)

  if (!is.null(brain_radius)) {
    brain_radius <- brain_radius * switch(trans_units, mm=1, cm=10, `in`=25.4)
  } else {
    brain_radius <- 50 # mm
  }

  rot_units <- match.arg(rot_units, rot_units)
  X[,4:6] <- X[,4:6] * switch(
    rot_units, 
    rad=brain_radius, deg=brain_radius*2*pi/360, 
    mm=1, cm=10, `in`=25.4
  )

  # Detrend if requested.
  if (detrend) {
    X <- dct_mat(X)
  }

  # Compute FD.
  Xdiff <- apply(X, 2, diff)
  FD <- c(0, apply(abs(Xdiff), 1, sum))

  # Revert units to `trans_units`.
  attr(FD, "units") <- trans_units
  switch(trans_units, mm=FD, cm=FD/10, `in`=FD/25.4)
}