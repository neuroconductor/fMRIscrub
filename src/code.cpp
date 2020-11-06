#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#include "tfdp.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


// [[Rcpp::export]]
List pcatf_core(arma::mat X, arma::mat U, arma::mat V, arma::vec d,
                double lambda, int K, int maxiter, int solveV,
                int verbose, double tol){

  int n = X.n_rows;
  int p = X.n_cols;

  arma::mat U_(n,K);
  arma::mat V_(p,K);
  arma::vec d_(K);


  arma::vec iters(K);
  arma::vec uk(n);
  arma::vec vk(p);
  arma::vec ukold(n);
  double dk;
  double err = 2*tol;


  for(int k = 0; k < K; k++){
    uk = U.col(k);
    dk = d(k);
    int iter=0;
    // vk = V.col(k);

    while(iter < maxiter){
      iter++;
      ukold = uk;
      Rcout << ukold(0) << std::endl;

      Rcout << n << std::endl;

      tf_dp(n, ukold, lambda);

      double ugh = arma::norm(ukold-uk);

      Rcout << ugh << std::endl;

      uk = arma::normalise(uk);

      vk = X.t()*uk;
      vk = arma::normalise(vk);
      err = arma::norm((uk-ukold)/n);
      if (err < tol) break;
    }

    iters(k) = iter;
    dk = arma::as_scalar( uk.t()*X*vk );
    X -= dk* uk * vk.t();
    U_.col(k) = uk;
    d_(k) = dk;
    if(solveV) V_.col(k) = vk;
  }

  List out = List::create(
    Named("U") = U_,
    Named("d") = d_,
    Named("V") = V_,
    Named("iters") = iters
  );
  return(out);
}
