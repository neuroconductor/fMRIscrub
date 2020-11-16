#include <RcppArmadillo.h>


// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::vec tf_dp(int n, arma::vec y, double lam) {
  int i;
  int k;
  int l;
  int r;
  int lo;
  int hi;
  double afirst;
  double alast;
  double bfirst;
  double blast;
  double alo;
  double blo;
  double ahi;
  double bhi;
  arma::vec x(2*n);
  arma::vec a(2*n);
  arma::vec b(2*n);
  arma::vec beta(n);

  /* These are the knots of the back-pointers */
  arma::vec tm(n-1);
  arma::vec tp(n-1);

  /* Take care of a few trivial cases */
  if(n==0) return y;
  if(n==1 || lam==0){
    for(i=0; i<n; i++) beta(i) = y(i);
    return beta;
  }

  // Rcout << "nontrivial at least" << std::endl;

  /* We step through the first iteration manually */
  tm(0) = -lam+y(0);
  tp(0) = lam+y(0);
  l = n-1;
  r = n;
  x(l) = tm(0);
  x(r) = tp(0);
  a(l) = 1;
  b(l) = -y(0)+lam;
  a(r) = -1;
  b(r) = y(0)+lam;
  afirst = 1;
  bfirst = -y(1)-lam;
  alast = -1;
  blast = y(1)-lam;

  /* Now iterations 2 through n-1 */
  for(k=1; k<n-1; k++){
    // Rcout << k << std::endl;
    /* Compute lo: step up from l until the derivative is greater than -lam */
    alo = afirst;
    blo = bfirst;
    for(lo=l; lo<=r; lo++){
      if(alo*x(lo)+blo > -lam) break;
      alo += a(lo);
      blo += b(lo);
    }

    /* Compute hi: step down from r until the derivative is less than lam */
    ahi = alast;
    bhi = blast;
    for(hi=r; hi>=lo; hi--) {
      if(-ahi*x(hi)-bhi < lam) break;
      ahi += a(hi);
      bhi += b(hi);
    }

    /* Compute the negative knot */
    tm(k) = (-lam-blo)/alo;
    l = lo-1;
    x(l) = tm(k);

    /* Compute the positive knot */
    tp(k) = (lam+bhi)/(-ahi);
    r = hi+1;
    x(r) = tp(k);

    /* Update a and b */
    a(l) = alo;
    b(l) = blo+lam;
    a(r) = ahi;
    b(r) = bhi+lam;
    afirst = 1;
    bfirst = -y(k+1)-lam;
    alast = -1;
    blast = y(k+1)-lam;
  }

  /* Compute the last coefficient: this is where the function has zero derivative */
  alo = afirst;
  blo = bfirst;
  for(lo=l; lo<=r; lo++){
    if(alo*x(lo)+blo > 0) break;
    alo += a(lo);
    blo += b(lo);
  }
  beta(n-1) = -blo/alo;

  /* Compute the rest of the coefficients, by the back-pointers */
  for(k=n-2; k>=0; k--){
    if(beta(k+1)>tp(k)) beta(k) = tp(k);
    else if(beta(k+1)<tm(k)) beta(k) = tm(k);
    else beta(k) = beta(k+1);
  }

  return beta;
}



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
  int iter;


  for(int k = 0; k < K; k++){
    uk = U.col(k);
    dk = d(k);
    vk = V.col(k);

    iter = 0;
    while(iter < maxiter){
      iter++;
      ukold = uk;
      uk = X*vk;
      uk = tf_dp(n, uk, lambda);
      uk = arma::normalise(uk);

      vk = X.t()*uk;
      vk = arma::normalise(vk);
      err = arma::norm((uk-ukold)/sqrt(n));
      if(err < tol) break;
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
