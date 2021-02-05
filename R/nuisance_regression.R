#' Nuisance regression
#'
#' Performs nuisance regression. The data and design matrix must both be
#'  centered, or an intercept must be included in the design matrix!
#'
#' @param X The TxV or VxT data.
#' @param design The TxQ matrix of nuisance regressors
#'
#' @return The data after nuisance regression
#' 
#' @export
nuisance_regression <- function(X, design){
  # https://stackoverflow.com/questions/19100600/extract-maximal-set-of-independent-columns-from-a-matrix
  # https://stackoverflow.com/questions/39167204/in-r-how-does-one-extract-the-hat-projection-influence-matrix-or-values-from-an
  qrd <- qr(design)
  design <- design[, qrd$pivot[seq_len(qrd$rank)]]
  qrd <- qr(design)
  Qd <- qr.Q(qrd)
  I_m_H <- diag(nrow(design)) - (Qd %*% t(Qd))
  if (nrow(X)==nrow(design)) {
    return(I_m_H %*% X)
  } else if (ncol(X)==nrow(design)) {
    return(X %*% I_m_H)
  } else {
    stop("X and design are not of compatible dimensions.")
  }
}