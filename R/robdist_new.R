#' Identify outliers based on bootstrap robust distance.
#' 
#' @param U N x P Dependent data
#' @param R_true N x N correlation matrix. If \code{NULL} (default), the identity
#'  matrix wil lbe used
#' @param boot_sample Number of bootstrap samples to compute. If \code{NULL},
#'  use the lowest reasonable value (\code{10/(1-quantile_cutoff)}).
#' @param quantile_cutoff The F-distribution quantile cutoff. Default: 
#'  \code{NA} (do not compute outliers).
#' 
#' @return List of "meas" (original RD values), "cut" (P cutoff values, one for 
#'  each observation) and "flag" (P indicators of outlyingness). If \code{is.null(quantile_cutoff)},
#'  the latter two entries are omitted.
#' 
#' @importFrom expm sqrtm
#' 
out_measures.robdist_bootstrap <- function(U, R_true=NULL, boot_sample=NULL, quantile_cutoff=NULL){
  if (is.null(R_true)) { R_true <- diag(nrow(U)) }

  x <- induc_indep(U, R_true)
  Y_tilde <- x$indep_data
  R_sqrt <- x$R_sqrt
  n <- nrow(Y_tilde)
  p <- ncol(Y_tilde)

  quant_cut2 <- quantile_cutoff # TO DO: take a look at this...
  if (is.null(quant_cut2)) { quant_cut2 <- 0.9999 }
  out1 <- out_measures.robdist(Y_tilde, quantile_cutoff=quant_cut2)
  best <- which(out1$info$inMCD)
  h <- length(best)
  
  # Obtain "notout"
  notbest <- setdiff(1:n, best)
  id_out <- which(out1$flag)
  notout <- setdiff(notbest, id_out)
  
  Y_tilde_in <- Y_tilde[best,]
  xbar_star <- colMeans(Y_tilde_in) # MCD estimate of mean , length of p
  Y_tilde_ins <- scale(Y_tilde_in, center = TRUE, scale = FALSE )
  nY_tilde <- nrow(Y_tilde_ins)
  S_star <- (t(Y_tilde_ins) %*% Y_tilde_ins)/(nY_tilde-1) # MCD estimate of covariance
  invcov <- solve(S_star) # dim p by p
  
  # Robust Distance of the Original Data
  RD_orig <- apply(U, 1, function(x) {t(x-xbar_star) %*% invcov %*% (x-xbar_star)})
  out <- list(meas=RD_orig)

  min_boot_sample <- ceiling(10 / (1 - quantile_cutoff))
  if (is.null(boot_sample)) { 
    boot_sample <- min_boot_sample 
  } else {
    if (boot_sample < min_boot_sample) {
      warning("The number of bootstrap samples should be at least `10/(1-quantile_cutoff)`. Using this value.\n")
      boot_sample <- min_boot_sample 
    }
  }

  if (!is.null(quantile_cutoff)){
    # pre-compute these to use for quickly computing robust distances for many observations simultaneously with matrix operations
    invcov_sqrt <- sqrtm(invcov)
    xbar_star_mat <- matrix(xbar_star, nrow=n, ncol=p , byrow = TRUE)
    
    # QUANTILE ADJUSTMENT for the "mostly notbest" observations
    gamma <- length(notbest)/ n
    alpha <- 1- quantile_cutoff
    adj_alpha <- alpha * gamma
    ### BOOTSTRAP STEP
    Y_tilde_boot <- array(NA, dim=c(n,p))
    RD_boot <- matrix(NA, n, boot_sample)
    for(b in 1:boot_sample){ # start loop over bootstrap samples
      Y_tilde_boot <- 0 * Y_tilde_boot
      best_boot <- sample(best, h, replace = T)
      notout_boot <- sample(notout, (n-h), replace=T)
      Y_tilde_boot[best,] <- Y_tilde[best_boot,]
      Y_tilde_boot[notbest,] <- Y_tilde[notout_boot,]
      Y_boot_b <- R_sqrt %*% Y_tilde_boot #RE-INDUCING DEPENDENCE
      #RD_boot[,b] <- apply(Y_boot_b,1, function(x) t(x-xbar_star) %*% invcov %*% (x-xbar_star))
      #following two lines are alternative way to get RDs
      temp <- (Y_boot_b - xbar_star_mat) %*% invcov_sqrt
      RD_boot[,b] <- rowSums(temp * temp)
    } # end loop over bootstrap samples
    RD_scaled <- out1$info$outMCD_scale * RD_boot # scale is from Hardin & Rocke's approach (2005)
    out$cut <- apply(RD_scaled[notbest,], 1, quantile, probs=c(1-adj_alpha))
    out$flag <- rep(FALSE, n)
    out$flag[notbest] <- (RD_orig[notbest] > out$cut)
  }

  out
}


#' Induce independence
#' 
#' @param dep_data N x P Dependent data
#' @param R_true N x N correlation matrix
#' 
#' @return List of "indep_data" (independent data) and "R_sqrt" (square root of R_true)
#' 
induc_indep <- function(dep_data, R_true){

  # TO DO
  # R_true
  #   1. Remove trends in the data

  R_svd <- svd(R_true)
  U <- R_svd$u
  V <- R_svd$v
  D_invsqrt <- diag(1/sqrt(R_svd$d))
  R_invsqrt <- U %*% D_invsqrt %*% t(V) # inverse square root of correlation matrix
  Y_tilde <- R_invsqrt %*% dep_data # n by p
  
  # for re-inducing dependence later
  D_sqrt <- diag(sqrt(R_svd$d))
  R_sqrt <- U %*% D_sqrt %*% t(V)
  
  list(indep_data=Y_tilde, R_sqrt=R_sqrt)
}