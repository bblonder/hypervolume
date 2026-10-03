#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix fastPdist2(NumericMatrix Ar, NumericMatrix Br) {
  int m = Ar.nrow(), 
  n = Br.nrow(),
  k = Ar.ncol();
  arma::mat A = arma::mat(Ar.begin(), m, k, false); 
  arma::mat B = arma::mat(Br.begin(), n, k, false); 
  
  arma::colvec An =  sum(square(A),1);
  arma::colvec Bn =  sum(square(B),1);
  
  arma::mat C = -2 * (A * B.t());
  C.each_col() += An;
  C.each_row() += Bn.t();

  // The |a|^2 + |b|^2 - 2 a.b expansion is fast but numerically unstable
  // for nearby points: rounding can push the squared distance slightly
  // negative and sqrt() of that returns NaN. hypervolume_distance(type =
  // "minimum") then hid the NaN with na.rm = TRUE, which silently changed
  // the answer for the exact case most at risk. Clamp negatives to zero
  // so the sqrt is always real. Values already >= 0 are untouched.
  C.transform([](double v) { return v < 0.0 ? 0.0 : v; });

  return wrap(sqrt(C));
}