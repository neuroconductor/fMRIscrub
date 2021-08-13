#' Summarize a \code{"scrub"} object
#'
#' Summary method for class \code{"scrub"}
#'
#' @param object Object of class \code{"scrub"}. 
#' @param ... further arguments passed to or from other methods.
#' @export
#' @method summary scrub
summary.scrub <- function(object, ...) {
  plot(object, ...)
}

#' @rdname summary.scrub
#' @export
#' 
#' @param x Object of class \code{"scrub"}. 
#' @method print summary.scrub
print.summary.scrub <- function(x, ...) {
  print(summary.scrub(x))
}

#' @rdname summary.scrub
#' @export
#' 
#' @method print scrub
print.scrub <- function(x, ...) {
  print.summary.scrub(summary(x, ...))
}