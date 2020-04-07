// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// rcpp_loop
NumericMatrix rcpp_loop(NumericMatrix mat_in, NumericMatrix transition, int n);
RcppExport SEXP _hermes6_rcpp_loop(SEXP mat_inSEXP, SEXP transitionSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type mat_in(mat_inSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type transition(transitionSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_loop(mat_in, transition, n));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_hermes6_rcpp_loop", (DL_FUNC) &_hermes6_rcpp_loop, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_hermes6(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
