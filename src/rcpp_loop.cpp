#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix rcpp_loop(NumericMatrix mat_in, NumericMatrix transition, int n) {
  NumericMatrix out(mat_in.nrow(), transition.ncol());
  out(0, _) = mat_in(0, _);
  NumericVector rm1, cm2;
  for (int i = 1; i < n; i++) {
    rm1 = out(i - 1,_);
    for (size_t j = 0; j < transition.ncol(); ++j) {
      cm2 = transition(_,j);
      out(i,j) = std::inner_product(rm1.begin(), rm1.end(), cm2.begin(), 0.);              
    }
  }
  return(out);
}
