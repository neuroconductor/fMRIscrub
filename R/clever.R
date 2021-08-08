#' Projection scrubbing
#' 
#' Projects the data onto directions of outlyingness and then calculates 
#'  leverage to identify outliers in high-dimensional data.
#' 
#' @section References:
#'  \itemize{
#'    \item{Mejia, A. F., Nebel, M. B., Eloyan, A., Caffo, B. & Lindquist, M. A. PCA leverage: outlier detection for high-dimensional functional magnetic resonance imaging data. Biostatistics 18, 521-536 (2017).}
#'    \item{Pham, D., McDonald, D., Ding, L., Nebel, M. B. & Mejia, A. Projection scrubbing: a more effective, data-driven fMRI denoising method. (2021).}
#'  }
#' 
#' @inheritSection clever_order_of_operations Order of operations
#' @inheritParams clever_Params
#' @param projection Choose one of the following: \code{"PCA"}, 
#'  \code{"PCATF"}, or \code{"ICA"}. The directions of outlyingness will be 
#'  selected from the high-kurtosis components among the top \eqn{k}, where
#'  \eqn{k} is the number of PCs selected by PESEL or that are above-average
#'  variance (see the \code{PESEL} argument). 
# @param R_true The \eqn{T \times T} correlation matrix, if known. Used for the bootstrap
#  robust distance measure.
#' @param PESEL Use \code{\link[pesel]{pesel}} to select the components? Default:
#'  \code{TRUE}. Otherwise, use the number of above-average variance PCs. 
#' @return A \code{"clever"} object, i.e. a list with components
#' \describe{
#'  \item{measure}{A numeric vector of leverage values.}
#'  \item{outlier_cutoff}{The numeric outlier cutoff value (\code{cutoff} times the median leverage).}
#'  \item{outlier_flag}{A logical vector where \code{TRUE} indicates suspected outlier presence.}
#'  \item{mask}{
#'    A length \eqn{P} numeric vector corresponding to the data locations in \code{X}. Each value indicates whether the location was masked:
#'    \describe{
#'      \item{0}{The data location was not masked out.}
#'      \item{-1}{The data location was masked out, because it had at least one \code{NA} or \code{NaN} value.}
#'      \item{-2}{The data location was masked out, because it was constant.}
#'    }
#'  }
#'  \item{PCA}{
#'    If PCA was used, this will be a list with components:
#'    \describe{
#'      \item{U}{The \eqn{T \times Q} PC score matrix.}
#'      \item{D}{The standard deviation of each PC.}
#'      \item{V}{The \eqn{P \times Q} PC directions matrix. Included only if \code{get_dirs}.}
#'      \item{highkurt}{The length \code{Q} logical vector indicating scores of high kurtosis.}
#'      \item{nPCs_PESEL}{The number of PCs selected by PESEL.}
#'      \item{nPCs_avgvar}{The number of above-average variance PCs.}
#'    }
#'    where \code{Q} is the number of PCs selected by PESEL or of above-average variance (or the greater of the two if both were used).
#'  }
#'  \item{PCATF}{
#'    If PCATF was used, this will be a list with components:
#'    \describe{
#'      \item{U}{The \eqn{T \times Q} PC score matrix.}
#'      \item{D}{The standard deviation of each PC.}
#'      \item{V}{The \eqn{P \times Q} PC directions matrix. Included only if \code{get_dirs}}
#'    }
#'  }
#'  \item{ICA}{
#'    If ICA was used, this will be a list with components:
#'    \describe{
#'      \item{S}{The \eqn{P \times Q} source signals matrix. Included only if \code{get_dirs}} 
#'      \item{M}{The \eqn{T \times Q} mixing matrix.}
#'      \item{highkurt}{The length \code{Q} logical vector indicating mixing scores of high kurtosis.}
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