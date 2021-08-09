#' Summarize a \code{"clever"} object
#'
#' Summary method for class \code{"clever"}
#'
#' @param object Object of class \code{"clever"}. 
#' @param ... further arguments passed to or from other methods.
#' @export
#' @method summary clever
summary.clever <- function(object, ...) {
  plot(object, ...)
}

#' @rdname summary.clever
#' @export
#' 
#' @param x Object of class \code{"clever"}. 
#' @method print summary.clever
print.summary.clever <- function(x, ...) {
  print(summary.clever(x))
}

#' @rdname summary.clever
#' @export
#' 
#' @method print clever
print.clever <- function(x, ...) {
  print.summary.clever(summary(x, ...))
}