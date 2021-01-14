#' Data-driven scrubbing with \code{clever}
#' 
#' Calculates data-driven scrubbing measures (leverage and DVARS) and identifies
#'  outliers in high-dimensional data.
#' 
#' @param X Wide numeric data matrix (\eqn{T observations \times V variables}, \eqn{T << V}).
#'  For example, if \code{X} represents an fMRI run, \eqn{T} should be the number
#'  of timepoints and \eqn{V} should be the number of brainordinate vertices/voxels.
#' 
#' @param measure Choose one of the following scrubbing measures:
#' 
#'  \describe{
#'    \item{\code{"leverage"}}{(Default) Leverage scrubbing, which is based on
#'      projecting the data onto directions thought to express outlier 
#'      information.}
#'    \item{\code{"DVARS"}}{traditional DVARS}
#'    \item{\code{"DVARS2"}}{Delta-percent-DVARS and z-score-DVARS (Afyouni and 
#'      Nichols, 2018)}
#'  }
#' 
#' @param projection Only applies if \code{measure} is \code{"leverage"}. Choose
#'  one of the following: \code{"PCA"}, \code{"PCATF"}, or \code{"ICA"}. The
#'  directions of outlyingness will be selected from components of this projection.
#' 
#' @param PESEL Only applies if \code{measure} is \code{"leverage"}. Leverage
#'  scrubbing selects directions among the largest \eqn{k} components. \eqn{k}
#'  can be determined by PESEL, a data-driven algorithm (default:
#'  \code{PESEL=TRUE}). Otherwise, \eqn{k} will be the number of
#'  above-average-variance PCs (\code{PESEL=FALSE}). 
#' 
#'  Both methods are based on the principal components (PCs), even if the
#'  projection is PCATF or ICA.
#' 
#'  For PCA and ICA, the final subset of components will be those in the top
#'  \eqn{k} which also have scores of high kurtosis (> 95th percentile of the
#'  kurtosis of data of equal length from a Normal distribution.)
#' @param solve_dirs Only applies if \code{measure} is \code{"leverage"}. 
#'  Should the projection directions be computed? Default: \code{FALSE}. 
#'  This will save memory, especially for PCA since the full SVD
#'  can be avoided. However, \code{solve_dirs=TRUE} is required to compute the
#'  leverage images.
#' @param center,scale Center the columns of the data by their medians, and scale the
#'  columns of the data by their median absolute distances (MADs)? Default: \code{TRUE}. 
#'  Centering is necessary for detrending and for computing PCA/ICA, so if this 
#'  is set to \code{FALSE}, the input data must already be centered.
#' @param DCT Detrend the columns of the data using the discrete cosine
#'  transform (DCT)? Use an integer to indicate the number of cosine bases to 
#'  use for detrending. Use \code{0} (default) to forgo detrending. 
#' 
#'  The data must be centered, either before input or with \code{center}.
#' 
#'  Detrending is highly recommended for time-series data, especially if there 
#'  are many time points or evolving circumstances affecting the data. Additionally,
#'  if kurtosis is being used to select the projection directions, trends can 
#'  induce positive or negative kurtosis, contaminating the connection between 
#'  high kurtosis and outlier presence. 
#'  
#'  Detrending should not be used with non-time-series data because the 
#'  observations are not temporally related.
#' @param nuisance_too A matrix of nuisance signals to regress from the data
#'  before, i.e. a "design matrix." Should have \eqn{T} rows. Nuisance
#'  regression will be performed simultaneously with DCT detrending if 
#'  applicable. \code{NULL} (default) to not add additional nuisance regressors.
#' @param PCATF_kwargs Arguments to \code{\link{PCATF}} in list form. Valid
#'  entries are:
#'  
#'  \describe{
#'    \item{K}{Maximum number of PCs to compute. Default: \code{100}.}
#'    \item{lambda}{Trend-filtering parameter. Default: \code{5}.}
#'    \item{niter_max}{Maximum number of iterations. Default: \code{1000}.}
#'    \item{verbose}{Print updates? Default: \code{FALSE}.}
#'  }
#' @param kurt_quantile Only applies to PCA and ICA leverage. What cutoff quantile
#'  for kurtosis should be used to select the components? Default: \code{0.95}.
#' @param get_outliers Should outliers be flagged based on cutoffs? Default: \code{TRUE}.
#' @param outlier_cutoff A number controlling the value of the scrubbing measure
#'  above which observations are labeled as outliers. If \code{NULL}, use a
#'  default:
#' 
#'  \describe{
#'    \item{\code{"leverage"}}{Minimum leverage value, in multiples of the median leverage. Default: \code{4} (will flag leverage scores more than four times the median).}
#'    \item{\code{"DVARS"}}{Minimum traditional-DVARS value. Default: \code{5}}
#'    \item{\code{"DVARS2"}}{A length-2 numeric vector representing a dual 
#'      Delta-percent-DVARS and z-score-DVARS cutoff. Both must be met for a 
#'      timepoint to be flagged. The Delta-percent-DVARS cutoff should be given
#'      in percentages; the z-score-DVARS cutoff should be given as a z-score. 
#'      Or, set this to "Afyouni-Nichols" (default) to require a
#'      Delta-percent-DVARS of more than 5\% and a z-score-DVARS greater than
#'      the right-tail 5\% significance level with Bonferroni FWER correction).}
#'  }
# @param R_true The \eqn{T \times T} correlation matrix, if known. Used for the bootstrap
#  robust distance measure.
#' @param full_PCA Return the full SVD? Default: \code{FALSE} (return
#'  only the components used to compute the measures).
#' @param verbose Should occasional updates be printed? Default: \code{FALSE}.
#'
#' @return A \code{"clever"} object, i.e. a list with components
#' \describe{
#'  \item{measure}{
#'    The scrubbing measure. If \code{"DVARS2"}, it will be a \code{data.frame}
#'    with the first column being DPDVARS and the second being ZDVARS.
#'  }
#'  \item{outlier_cutoff}{
#'    The value of the scrubbing measure used to identify outliers (i.e.
#'    minimum leverage, traditional DVARS, or DPDVARS and ZDVARS).
#'  }
#'  \item{outlier_flags}{
#'    A length \eqn{T} logical vector. \code{TRUE} values indicate suspected outliers. 
#'  }
#'  \item{PCA}{
#'    If PCA was used, this will be a list with components:
#'    \describe{
#'      \item{U}{The \eqn{T \times Q} PC score matrix.}
#'      \item{D}{The standard deviation of each PC.}
#'      \item{V}{The \eqn{P \times Q} PC directions matrix. Included only if \code{solve_dirs}}
#'      \item{highkurt}{The length \code{Q} logical vector indicating scores of high kurtosis.}

#'    }
#'  }
#'  \item{PCATF}{
#'    If PCATF was used, this will be a list with components:
#'    \describe{
#'      \item{U}{The \eqn{T \times Q} PC score matrix.}
#'      \item{D}{The standard deviation of each PC.}
#'      \item{V}{The \eqn{P \times Q} PC directions matrix. Included only if \code{solve_dirs}}
#'    }
#'  }
#'  \item{ICA}{
#'    If ICA was used, this will be a list with components:
#'    \describe{
#'      \item{S}{The \eqn{P \times Q} source signals matrix.} 
#'      \item{M}{The \eqn{T \times Q} mixing matrix.}
#'      \item{highkurt}{The length \code{Q} logical vector indicating scores of high kurtosis.}
#'    }
#'  }
#'  \item{call}{The call to this function.}
#' }
#'
#' @importFrom pesel pesel
#' @importFrom robustbase rowMedians
#' @importFrom stats mad qnorm var setNames
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
  X, measure="leverage", projection="PCA", PESEL=TRUE, solve_dirs=FALSE,
  center=TRUE, scale=TRUE, DCT=0, nuisance_too=NULL,
  PCATF_kwargs=NULL, kurt_quantile=.95,
  get_outliers=TRUE, 
  outlier_cutoff=NULL,
  full_PCA=FALSE,
  verbose=FALSE){
  
  # Check and format any arguments handled differently than `clever_multi`.
  measure <- match.arg(measure, c("leverage", "DVARS", "DVARS2"))
  stopifnot(length(measure) == 1)
  stopifnot(is.logical(PESEL))
  projection <- match.arg(projection, c("PCA", "ICA", "PCATF"))
  projection <- switch(projection,
    PCA = ifelse(PESEL, "PCA2_kurt", "PCA_kurt"),
    ICA = ifelse(PESEL, "ICA2_kurt", "ICA_kurt"),
    PCATF = "PCATF"
  )

  if (!is.null(outlier_cutoff)) {
    outlier_cutoff <- list(outlier_cutoff)
    names(outlier_cutoff) <- measure
  } else {
    outlier_cutoff <- switch(measure,
      leverage = 4,
      DVARS  = 5,
      DVARS2 = "Afyouni-Nichols"
    )
  }

  # Run `clever_multi`.
  clev <- clever_multi(
    X=X, measures=measure, projections=projection, solve_dirs=solve_dirs,
    center=center, scale=scale, DCT=DCT, nuisance_too=nuisance_too, 
    PCATF_kwargs=PCATF_kwargs, kurt_quantile=kurt_quantile, 
    get_outliers=get_outliers, outlier_cutoffs=outlier_cutoff,
    full_PCA=full_PCA, verbose=verbose
  )

  # Re-format results.
  clever_from_multi(clev)
}