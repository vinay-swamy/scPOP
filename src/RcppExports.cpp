// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// compute_simpson_index
arma::vec compute_simpson_index(arma::mat& D, arma::umat& knn_idx, arma::vec& batch_labels, int n_batches, double perplexity, double tol);
RcppExport SEXP _scPOP_compute_simpson_index(SEXP DSEXP, SEXP knn_idxSEXP, SEXP batch_labelsSEXP, SEXP n_batchesSEXP, SEXP perplexitySEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type D(DSEXP);
    Rcpp::traits::input_parameter< arma::umat& >::type knn_idx(knn_idxSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type batch_labels(batch_labelsSEXP);
    Rcpp::traits::input_parameter< int >::type n_batches(n_batchesSEXP);
    Rcpp::traits::input_parameter< double >::type perplexity(perplexitySEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_simpson_index(D, knn_idx, batch_labels, n_batches, perplexity, tol));
    return rcpp_result_gen;
END_RCPP
}
// countPairs
List countPairs(IntegerVector classi1, IntegerVector classi2, IntegerVector order);
RcppExport SEXP _scPOP_countPairs(SEXP classi1SEXP, SEXP classi2SEXP, SEXP orderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type classi1(classi1SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type classi2(classi2SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type order(orderSEXP);
    rcpp_result_gen = Rcpp::wrap(countPairs(classi1, classi2, order));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_scPOP_compute_simpson_index", (DL_FUNC) &_scPOP_compute_simpson_index, 6},
    {"_scPOP_countPairs", (DL_FUNC) &_scPOP_countPairs, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_scPOP(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
