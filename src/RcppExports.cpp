// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// DFS
Rcpp::IntegerVector DFS(SEXP amattt, SEXP n_nodes);
RcppExport SEXP _MRSL_DFS(SEXP amatttSEXP, SEXP n_nodesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type amattt(amatttSEXP);
    Rcpp::traits::input_parameter< SEXP >::type n_nodes(n_nodesSEXP);
    rcpp_result_gen = Rcpp::wrap(DFS(amattt, n_nodes));
    return rcpp_result_gen;
END_RCPP
}
// STEP1
SEXP STEP1(SEXP data_sum_beta, SEXP data_sum_se, SEXP beta1, SEXP use_eggers_step1);
RcppExport SEXP _MRSL_STEP1(SEXP data_sum_betaSEXP, SEXP data_sum_seSEXP, SEXP beta1SEXP, SEXP use_eggers_step1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type data_sum_beta(data_sum_betaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type data_sum_se(data_sum_seSEXP);
    Rcpp::traits::input_parameter< SEXP >::type beta1(beta1SEXP);
    Rcpp::traits::input_parameter< SEXP >::type use_eggers_step1(use_eggers_step1SEXP);
    rcpp_result_gen = Rcpp::wrap(STEP1(data_sum_beta, data_sum_se, beta1, use_eggers_step1));
    return rcpp_result_gen;
END_RCPP
}
// STEP2
SEXP STEP2(SEXP amattt, SEXP data_sum_beta, SEXP data_sum_se, SEXP beta1, SEXP adj_methods, SEXP use_eggers_step2, SEXP vary_mvmr_adj);
RcppExport SEXP _MRSL_STEP2(SEXP amatttSEXP, SEXP data_sum_betaSEXP, SEXP data_sum_seSEXP, SEXP beta1SEXP, SEXP adj_methodsSEXP, SEXP use_eggers_step2SEXP, SEXP vary_mvmr_adjSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type amattt(amatttSEXP);
    Rcpp::traits::input_parameter< SEXP >::type data_sum_beta(data_sum_betaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type data_sum_se(data_sum_seSEXP);
    Rcpp::traits::input_parameter< SEXP >::type beta1(beta1SEXP);
    Rcpp::traits::input_parameter< SEXP >::type adj_methods(adj_methodsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type use_eggers_step2(use_eggers_step2SEXP);
    Rcpp::traits::input_parameter< SEXP >::type vary_mvmr_adj(vary_mvmr_adjSEXP);
    rcpp_result_gen = Rcpp::wrap(STEP2(amattt, data_sum_beta, data_sum_se, beta1, adj_methods, use_eggers_step2, vary_mvmr_adj));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MRSL_DFS", (DL_FUNC) &_MRSL_DFS, 2},
    {"_MRSL_STEP1", (DL_FUNC) &_MRSL_STEP1, 4},
    {"_MRSL_STEP2", (DL_FUNC) &_MRSL_STEP2, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_MRSL(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}