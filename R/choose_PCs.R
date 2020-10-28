#' Identifies the components of sufficient kurtosis.
#'
#' Each component is detrended (this can be disabled) before computing kurtosis. 
#'  The kurtosis cutoff is the 90% quantile of the sampling distribution of 
#'  kurtosis for Normal data of the same length as the components; it is
#'  estimated by simulation or calculated from the theoretical asymptotic
#'  distribution if the components are long enough.
#'
#' @param Comps A matrix; each column is a component. For PCA, this is the U
#'  matrix. For ICA, this is the M matrix.
#' @param kurt_quantile components with kurtosis of at least this quantile are kept.
#' @param detrend Should components be detrended before measuring kurtosis? Default is
#'  \code{TRUE}. Recommended if observations represent a time series. Use
#'  \code{FALSE} if they have already been detrended.
#' @param n_sim The number of simulation data to use for estimating the sampling
#'  distribution of kurtosis. Only used if a new simulation is performed. (If
#'  \eqn{n<1000} and the quantile is 90%, a pre-computed value is used instead.
#'  If \eqn{n>1000}, the theoretical asymptotic distribution is used instead.)
#'
#' @return A logical vector indicating whether each component has high kurtosis.
#'
#' @importFrom stats quantile qnorm
#' @importFrom e1071 kurtosis
#' @importFrom MASS mvrnorm
#' @export
high_kurtosis <- function(Comps, 
  kurt_quantile = 0.9, detrend = TRUE, n_sim = 5000){

  m <- nrow(Comps); n <- ncol(Comps)

  if (detrend) {  Comps <- Comps - apply(Comps, 2, est_trend) }

  kurt <- apply(Comps, 2, kurtosis, type=1)

  # Determine the quantile cutoff.
  if(m < 1000){
    if(kurt_quantile == .9){
      # Use precomputed empirical 0.90 quantile.
      cut <- kurt_90_quant[m]
    } else {
      # Simulate and compute the empirical quantile if not 0.90.
      sim <- apply(t(mvrnorm(n_sim, mu=rep(0, m), diag(m))), 2, kurtosis, type=1)
      cut <- quantile(sim, kurt_quantile)
    }
  } else {
    # Use theoretical quantile.
    cut <- qnorm(kurt_quantile) * sqrt( (24*m*(m-1)^2) / ((m-3)*(m-2)*(m+3)*(m+5)) )
  }

  # Constant components will have NaN kurtosis.
  kurt[kurt %in% c(NA, NaN)] <- -Inf

  kurt > cut
}
