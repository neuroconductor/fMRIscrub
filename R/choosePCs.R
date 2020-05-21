#' Selects the principle components of above-average variance from a SVD.
#'
#' PCs with above-average variance are retained, but the total number kept is
#'  constrained by the \code{max_keep} and \code{min_keep} arguments.
#'
#' @param svd An SVD decomposition; i.e. a list containing u, d, and v.
#' @param max_keep If specified, the total number kept will be at most this
#'  value.
#' @param min_keep The total number kept will be at least this
#'  value. The default value is \code{1}.
#'
#' @return The original indices of the PCs which were retained, in order of
#'  decreasing variance (i.e. increasing index).
#'
#' @export
choose_PCs.variance <- function(svd, max_keep = NULL, min_keep = 1){
  var <- svd$d^2

  # Identify how many PCs will be kept.
  n_keep <- sum(var > mean(var))

  # Constrain this number kept between the minumum and maximum, if specified.
  if(!is.null(max_keep)){n_keep <- min(n_keep, max_keep)}
  if(!is.null(min_keep)){n_keep <- max(n_keep, min_keep)}

  # PCs are already ordered by decreasing variance.
  indices <- 1:n_keep

  return(indices)
}

#' Selects the principle components (PCs) of sufficient kurtosis from a SVD.
#'
#' First, the largest PCs which together explain 90% of the variance are
#'  retained, and smaller components are removed. Each PC is detrended (this can
#'  be disabled). The kurtosis cutoff is then the 90% quantile of the sampling
#'  distribution of kurtosis for Normal data of the same length as the PCs; it
#'  is estimated by simulation or calculated from the theoretical asymptotic
#'  distribution if the PCs are long enough. Finally, the total number kept is
#'  constrained by the \code{max_keep} and \code{min_keep} arguments.
#'
#' @param svd An SVD decomposition; i.e. a list containing u, d, and v.
#' @param kurt_quantile PCs with kurtosis of at least this quantile are kept.
#' @param detrend Should PCs be detrended before measuring kurtosis? Default is
#'  \code{TRUE}. Recommended if observations represent a time series.
#' @param max_keep If specified, the total number kept will be at most this
#'  value.
#' @param min_keep The total number kept will be at least this
#'  value. The default value is \code{1}.
#' @param n_sim The number of simulation data to use for estimating the sampling
#'  distribution of kurtosis. Only used if a new simulation is performed. (If
#'  \eqn{n<1000} and the quantile is 90%, a pre-computed value is used instead.
#'  If \eqn{n>1000}, the theoretical asymptotic distribution is used instead.)
#'
#' @return The original indices of the PCs which were retained, in order of
#'  decreasing kurtosis.
#'
#' @importFrom e1071 kurtosis
#' @importFrom MASS mvrnorm
#' @export
choose_PCs.kurtosis <- function(svd, kurt_quantile = 0.9, detrend = TRUE,
  max_keep = NULL, min_keep = 1, n_sim = 5000){
  U <- svd$u
  m <- nrow(U)

  # First remove components that explain less than 90% of variation.
  cumvarexp <- cumsum(svd$d/sum(svd$d))
  n <- max(min(which((cumvarexp > .90))), min_keep)
  n <- max(n, min_keep)
  U <- U[,1:n]
  if(n==1){U <- matrix(U, ncol=1)}

  # Compute the kurtosis of the remaining PCs, detrending if applicable.
  if(detrend){
    U.dt <- U - apply(U, 2, est_trend)
    kurt <- apply(U.dt, 2, kurtosis, type=1)
  } else {
    kurt <- apply(U, 2, kurtosis, type=1)
  }

  # Determine the quantile cutoff.
  if(m < 1000){
    if(kurt_quantile == .9){
      # Use precomputed empirical quantile.
      cut <- kurt_90_quant[m]
    } else {
      # Simulate and compute the empirical quantile.
      sim <- apply(t(mvrnorm(n_sim, mu=rep(0, m), diag(m))), 2, kurtosis, type=1)
      cut <- quantile(sim, kurt_quantile)
    }
  } else {
    # Use theoretical quantile.
    cut <- qnorm(kurt_quantile) * sqrt( (24*m*(m-1)^2) / ((m-3)*(m-2)*(m+3)*(m+5)) )
  }

  # Identify how many PCs will be kept.
  n_keep <- sum(kurt > cut)

  # Constrain the number kept between the minumum and maximum, if specified.
  if(!is.null(max_keep)){n_keep <- min(n_keep, max_keep)}
  if(!is.null(min_keep)){n_keep <- max(n_keep, min_keep)}

  # The PCs with greatest kurtosis are chosen.
  indices <- order(-kurt)[1:n_keep]

  return(indices)
}
