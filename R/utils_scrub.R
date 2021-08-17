#' Summarize a \code{"scrub_projection"} object
#'
#' Summary method for class \code{"scrub_projection"}
#'
#' @param object Object of class \code{"scrub_projection"}. 
#' @param ... further arguments passed to or from other methods.
#' @return A plot of the scrubbing results
#' @export
#' @method summary scrub_projection
summary.scrub_projection <- function(object, ...) {
  plot(object, ...)
}

#' @rdname summary.scrub_projection
#' @export
#' 
#' @param x Object of class \code{"scrub_projection"}. 
#' @method print summary.scrub_projection
print.summary.scrub_projection <- function(x, ...) {
  print(summary.scrub_projection(x))
}

#' @rdname summary.scrub_projection
#' @export
#' 
#' @method print scrub_projection
print.scrub_projection <- function(x, ...) {
  print.summary.scrub_projection(summary(x, ...))
}

#' Summarize a \code{"scrub_DVARS"} object
#'
#' Summary method for class \code{"scrub_DVARS"}
#'
#' @param object Object of class \code{"scrub_DVARS"}. 
#' @param ... further arguments passed to or from other methods.
#' @export
#' @return A plot of the scrubbing results
#' @method summary scrub_DVARS
summary.scrub_DVARS <- function(object, ...) {
  plot(object, ...)
}

#' @rdname summary.scrub_DVARS
#' @export
#' 
#' @param x Object of class \code{"scrub_DVARS"}. 
#' @method print summary.scrub_DVARS
print.summary.scrub_DVARS <- function(x, ...) {
  print(summary.scrub_DVARS(x))
}

#' @rdname summary.scrub_DVARS
#' @export
#' 
#' @method print scrub_DVARS
print.scrub_DVARS <- function(x, ...) {
  print.summary.scrub_DVARS(summary(x, ...))
}

#' Summarize a \code{"scrub_FD"} object
#'
#' Summary method for class \code{"scrub_FD"}
#'
#' @param object Object of class \code{"scrub_FD"}. 
#' @param ... further arguments passed to or from other methods.
#' @export
#' @return A plot of the scrubbing results
#' @method summary scrub_FD
summary.scrub_FD <- function(object, ...) {
  plot(object, ...)
}

#' @rdname summary.scrub_FD
#' @export
#' 
#' @param x Object of class \code{"scrub_FD"}. 
#' @method print summary.scrub_FD
print.summary.scrub_FD <- function(x, ...) {
  print(summary.scrub_FD(x))
}

#' @rdname summary.scrub_FD
#' @export
#' 
#' @method print scrub_FD
print.scrub_FD <- function(x, ...) {
  print.summary.scrub_FD(summary(x, ...))
}