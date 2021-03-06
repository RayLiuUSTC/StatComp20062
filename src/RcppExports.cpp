// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// dLaplace
double dLaplace(double x);
RcppExport SEXP _StatComp20062_dLaplace(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(dLaplace(x));
    return rcpp_result_gen;
END_RCPP
}
// rw_Metropolis
List rw_Metropolis(double sigma, double x0, int N);
RcppExport SEXP _StatComp20062_rw_Metropolis(SEXP sigmaSEXP, SEXP x0SEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(rw_Metropolis(sigma, x0, N));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_StatComp20062_dLaplace", (DL_FUNC) &_StatComp20062_dLaplace, 1},
    {"_StatComp20062_rw_Metropolis", (DL_FUNC) &_StatComp20062_rw_Metropolis, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_StatComp20062(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
