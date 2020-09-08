#' Determine the leverage outliers
#' 
#' Identify outliers based on leverage values above a multiple of the median.
#' 
#' @param lev The leverage values
#' @param median_cutoff The cutoff, in multiples of the median of \code{lev}.
#'  Default: \code{4}
#' @return List with entries \code{"cut"} (the leverage cutoff value) and
#'  \code{"flag"} (logical vector indicating the outliers)
outs.leverage <- function(lev, median_cutoff=4){
  cut <- median_cutoff * median(lev)
  list(cut = cut, flag= lev > cut)
}

#' Determine the robust distance outliers
#' 
#' Identify outliers based on leverage values above a multiple of the median.
#' 
#' @param rbd The robust distance values
#' @param Fparam_1,Fparam_2 the F distribution degrees of freedom
#' @param inMCD Logical vector indicating in-MCD observations
#' @param outMCD_scale The scale for out-of-MCD observations
#' @param quantile_cutoff The F-distribution quantile cutoff. Default: 
#'  \code{.9999}
#' @return List with entries \code{"cut"} (the robust distance cutoff value for
#'  scaled out-of-MCD observations) and \code{"flag"} (logical vector 
#'  indicating the outliers)
outs.robdist <- function(
  rbd, Fparam_1, Fparam_2, inMCD, outMCD_scale, quantile_cutoff=.9999){
  cut <- qf(p=quantile_cutoff, df1=Fparam_1, df2=Fparam_2)
  list(cut = cut, flag = ifelse(inMCD, FALSE, rbd * outMCD_scale > cut))
}