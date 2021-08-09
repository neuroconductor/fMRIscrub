#' Projection scrubbing
#' 
#' Projection scrubbing is a data-driven method for identifying artifactual
#'  volumes in fMRI scans. It works by identifying component directions in the 
#'  data likely to represent artifactual patterns, and then computing a 
#'  composite measure for the amount of these signals at each volume. The 
#'  projection can be PCA, ICA, or "fused PCA," and the composite measure is 
#'  based on the linear model concept of leverage. Projection scrubbing can 
#'  also be used for other outlier detection tasks involving high-dimensional 
#'  data.
#' 
#'  Refer to the projection scrubbing vignette for a demonstration and an
#'  outline of the algorithm: \code{vignette("projection_scrubbing", package="cubature")}
#' 
#' @section References:
#'  \itemize{
#'    \item{Mejia, A. F., Nebel, M. B., Eloyan, A., Caffo, B. & Lindquist, M. A. PCA leverage: outlier detection for high-dimensional functional magnetic resonance imaging data. Biostatistics 18, 521-536 (2017).}
#'    \item{Pham, D., McDonald, D., Ding, L., Nebel, M. B. & Mejia, A. Projection scrubbing: a more effective, data-driven fMRI denoising method. (2021).}
#'  }
#' 
#' @inheritParams clever_Params
#' @param projection One of the following: \code{"ICA"} (default), \code{"PCA"}, 
#'  or \code{"PCATF"}.
# @param R_true The \eqn{T} by \eqn{T} correlation matrix, if known. Used for the bootstrap
#  robust distance measure.
#' @param PESEL Use \code{\link[pesel]{pesel}} to select the number of 
#'  components? Default: \code{TRUE}. Otherwise, use the number of principal
#'  components with above-average variance.
#' @return A \code{"clever"} object, i.e. a list with components
#' \describe{
#'  \item{measure}{A numeric vector of leverage values.}
#'  \item{outlier_cutoff}{The numeric outlier cutoff value (\code{cutoff} times the median leverage).}
#'  \item{outlier_flag}{A logical vector where \code{TRUE} indicates where leverage exceeds the cutoff, signalling suspected outlier presence.}
#'  \item{mask}{
#'    A length \eqn{P} numeric vector corresponding to the data locations in \code{X}. Each value indicates whether the location was masked:
#'    \describe{
#'      \item{0}{The data location was not masked out.}
#'      \item{-1}{The data location was masked out, because it had at least one \code{NA} or \code{NaN} value.}
#'      \item{-2}{The data location was masked out, because it was constant.}
#'    }
#'  }
#'  \item{PCA}{
#'    This will be a list with components:
#'    \describe{
#'      \item{U}{The \eqn{T} by \eqn{Q} PC score matrix.}
#'      \item{D}{The standard deviation of each PC.}
#'      \item{V}{The \eqn{P} by \eqn{Q} PC directions matrix. Included only if \code{get_dirs}.}
#'      \item{highkurt}{The length \code{Q} logical vector indicating scores of high kurtosis.}
#'      \item{U_dt}{Detrended components of \code{U}. Included only if components were mean- or variance-detrended.}
#'      \item{highkurt}{The length \code{Q} logical vector indicating detrended scores of high kurtosis.}
#'      \item{nPCs_PESEL}{The number of PCs selected by PESEL.}
#'      \item{nPCs_avgvar}{The number of above-average variance PCs.}
#'    }
#'    where \code{Q} is the number of PCs selected by PESEL or of above-average variance (or the greater of the two if both were used). 
#'    If PCA was not used, all entries except \code{nPCs_PESEL} and/or \code{nPCs_avgvar} will not be included, depending on which
#'    method(s) was used to select the number of PCs.
#'  }
#'  \item{PCATF}{
#'    If PCATF was used, this will be a list with components:
#'    \describe{
#'      \item{U}{The \eqn{T} by \eqn{Q} PC score matrix.}
#'      \item{D}{The standard deviation of each PC.}
#'      \item{V}{The \eqn{P} by \eqn{Q} PC directions matrix. Included only if \code{get_dirs}}
#'      \item{highkurt}{The length \code{Q} logical vector indicating scores of high kurtosis.}
#'      \item{U_dt}{Detrended components of \code{U}. Included only if components were mean- or variance-detrended.}
#'      \item{highkurt}{The length \code{Q} logical vector indicating detrended scores of high kurtosis. Included only if components were mean- or variance-detrended.}
#'    }
#'  }
#'  \item{ICA}{
#'    If ICA was used, this will be a list with components:
#'    \describe{
#'      \item{S}{The \eqn{P} by \eqn{Q} source signals matrix. Included only if \code{get_dirs}} 
#'      \item{M}{The \eqn{T} by \eqn{Q} mixing matrix.}
#'      \item{highkurt}{The length \code{Q} logical vector indicating mixing scores of high kurtosis.}
#'      \item{M_dt}{Detrended components of \code{M}. Included only if components were mean- or variance-detrended.}
#'      \item{highkurt}{The length \code{Q} logical vector indicating detrended mixing scores of high kurtosis. Included only if components were mean- or variance-detrended.}
#'    }
#'  }
#' }
#'
#' @export
#'
#' @examples
#' n_voxels = 1e4
#' n_timepoints = 100
#' X = matrix(rnorm(n_timepoints*n_voxels), ncol = n_voxels)
#'
#' clev = clever(X)
clever = function(
  X, projection=c("PCA", "PCATF", "ICA"), 
  nuisance=cbind(1, dct_bases(nrow(X), 4)),
  center=TRUE, scale=TRUE, comps_mean_dt=FALSE, comps_var_dt=FALSE,
  PESEL=TRUE, kurt_quantile=.99, PCATF_kwargs=NULL, 
  get_dirs=FALSE, full_PCA=FALSE,
  get_outliers=TRUE, cutoff=4, 
  verbose=FALSE){
  
  projection <- match.arg(projection, c("PCA", "PCATF", "ICA"))
  if (!PESEL) { projection <- paste0(projection, "2") }
  if (kurt_quantile > 0) { 
    projection <- paste0(projection, "_kurt")
  } else {
    kurt_quantile <- .99
  }

  # Run `clever_multi`.
  clev <- clever_multi(
    X=X, projection=projection, 
    nuisance=nuisance,
    center=center, scale=scale, comps_mean_dt=comps_mean_dt, comps_var_dt=comps_var_dt,
    kurt_quantile=kurt_quantile, PCATF_kwargs=PCATF_kwargs,
    get_dirs=get_dirs, full_PCA=full_PCA,
    get_outliers=get_outliers, cutoff=cutoff,
    verbose=verbose
  )

  # Re-format results.
  clever_from_multi(clev)
}