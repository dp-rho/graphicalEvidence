// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// print_times
void print_times(const int nruns);
RcppExport SEXP _graphicalEvidence_print_times(SEXP nrunsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type nruns(nrunsSEXP);
    print_times(nruns);
    return R_NilValue;
END_RCPP
}
// reset_times
void reset_times();
RcppExport SEXP _graphicalEvidence_reset_times() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    reset_times();
    return R_NilValue;
END_RCPP
}
// rgamma_compiled
NumericVector rgamma_compiled(int n, NumericVector vec_shapes, NumericVector vec_rates);
RcppExport SEXP _graphicalEvidence_rgamma_compiled(SEXP nSEXP, SEXP vec_shapesSEXP, SEXP vec_ratesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type vec_shapes(vec_shapesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type vec_rates(vec_ratesSEXP);
    rcpp_result_gen = Rcpp::wrap(rgamma_compiled(n, vec_shapes, vec_rates));
    return rcpp_result_gen;
END_RCPP
}
// bind_random_samples
void bind_random_samples(NumericVector rgammas, NumericVector rnorms);
RcppExport SEXP _graphicalEvidence_bind_random_samples(SEXP rgammasSEXP, SEXP rnormsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type rgammas(rgammasSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rnorms(rnormsSEXP);
    bind_random_samples(rgammas, rnorms);
    return R_NilValue;
END_RCPP
}
// bind_random_samples_rmatrix
void bind_random_samples_rmatrix(NumericVector rgammas, NumericVector rnorms, NumericVector rgigs, NumericVector runis);
RcppExport SEXP _graphicalEvidence_bind_random_samples_rmatrix(SEXP rgammasSEXP, SEXP rnormsSEXP, SEXP rgigsSEXP, SEXP runisSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type rgammas(rgammasSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rnorms(rnormsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rgigs(rgigsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type runis(runisSEXP);
    bind_random_samples_rmatrix(rgammas, rnorms, rgigs, runis);
    return R_NilValue;
END_RCPP
}
// mcmc_hw
List mcmc_hw(int n, int burnin, int nmc, int alpha, int p, NumericVector s_mat_nvec, NumericVector scale_mat_nvec, NumericVector g_mat_adj_nvec, NumericVector gibbs_mat_nvec, NumericVector init_gibbs_nvec);
RcppExport SEXP _graphicalEvidence_mcmc_hw(SEXP nSEXP, SEXP burninSEXP, SEXP nmcSEXP, SEXP alphaSEXP, SEXP pSEXP, SEXP s_mat_nvecSEXP, SEXP scale_mat_nvecSEXP, SEXP g_mat_adj_nvecSEXP, SEXP gibbs_mat_nvecSEXP, SEXP init_gibbs_nvecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< int >::type nmc(nmcSEXP);
    Rcpp::traits::input_parameter< int >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s_mat_nvec(s_mat_nvecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type scale_mat_nvec(scale_mat_nvecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type g_mat_adj_nvec(g_mat_adj_nvecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type gibbs_mat_nvec(gibbs_mat_nvecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type init_gibbs_nvec(init_gibbs_nvecSEXP);
    rcpp_result_gen = Rcpp::wrap(mcmc_hw(n, burnin, nmc, alpha, p, s_mat_nvec, scale_mat_nvec, g_mat_adj_nvec, gibbs_mat_nvec, init_gibbs_nvec));
    return rcpp_result_gen;
END_RCPP
}
// mcmc_hw_rmatrix
List mcmc_hw_rmatrix(int n, int burnin, int nmc, int p, int prior, int dof, double lambda, NumericVector s_mat_nvec, NumericVector gibbs_mat_nvec);
RcppExport SEXP _graphicalEvidence_mcmc_hw_rmatrix(SEXP nSEXP, SEXP burninSEXP, SEXP nmcSEXP, SEXP pSEXP, SEXP priorSEXP, SEXP dofSEXP, SEXP lambdaSEXP, SEXP s_mat_nvecSEXP, SEXP gibbs_mat_nvecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< int >::type nmc(nmcSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< int >::type dof(dofSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s_mat_nvec(s_mat_nvecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type gibbs_mat_nvec(gibbs_mat_nvecSEXP);
    rcpp_result_gen = Rcpp::wrap(mcmc_hw_rmatrix(n, burnin, nmc, p, prior, dof, lambda, s_mat_nvec, gibbs_mat_nvec));
    return rcpp_result_gen;
END_RCPP
}
// mcmc_last_col
List mcmc_last_col(const unsigned int n, const unsigned int burnin, const unsigned int nmc, const unsigned int alpha, const unsigned int p, NumericVector fixed_last_col_nvec, NumericVector s_mat_nvec, NumericVector scale_mat_nvec, NumericVector g_mat_adj_nvec, NumericVector gibbs_mat_nvec, NumericVector post_mean_omega_nvec);
RcppExport SEXP _graphicalEvidence_mcmc_last_col(SEXP nSEXP, SEXP burninSEXP, SEXP nmcSEXP, SEXP alphaSEXP, SEXP pSEXP, SEXP fixed_last_col_nvecSEXP, SEXP s_mat_nvecSEXP, SEXP scale_mat_nvecSEXP, SEXP g_mat_adj_nvecSEXP, SEXP gibbs_mat_nvecSEXP, SEXP post_mean_omega_nvecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const unsigned int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type nmc(nmcSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type p(pSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type fixed_last_col_nvec(fixed_last_col_nvecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s_mat_nvec(s_mat_nvecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type scale_mat_nvec(scale_mat_nvecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type g_mat_adj_nvec(g_mat_adj_nvecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type gibbs_mat_nvec(gibbs_mat_nvecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type post_mean_omega_nvec(post_mean_omega_nvecSEXP);
    rcpp_result_gen = Rcpp::wrap(mcmc_last_col(n, burnin, nmc, alpha, p, fixed_last_col_nvec, s_mat_nvec, scale_mat_nvec, g_mat_adj_nvec, gibbs_mat_nvec, post_mean_omega_nvec));
    return rcpp_result_gen;
END_RCPP
}
// set_cores
void set_cores(const int cores);
RcppExport SEXP _graphicalEvidence_set_cores(SEXP coresSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type cores(coresSEXP);
    set_cores(cores);
    return R_NilValue;
END_RCPP
}
// set_seed
void set_seed(unsigned int seed);
RcppExport SEXP _graphicalEvidence_set_seed(SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type seed(seedSEXP);
    set_seed(seed);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_graphicalEvidence_print_times", (DL_FUNC) &_graphicalEvidence_print_times, 1},
    {"_graphicalEvidence_reset_times", (DL_FUNC) &_graphicalEvidence_reset_times, 0},
    {"_graphicalEvidence_rgamma_compiled", (DL_FUNC) &_graphicalEvidence_rgamma_compiled, 3},
    {"_graphicalEvidence_bind_random_samples", (DL_FUNC) &_graphicalEvidence_bind_random_samples, 2},
    {"_graphicalEvidence_bind_random_samples_rmatrix", (DL_FUNC) &_graphicalEvidence_bind_random_samples_rmatrix, 4},
    {"_graphicalEvidence_mcmc_hw", (DL_FUNC) &_graphicalEvidence_mcmc_hw, 10},
    {"_graphicalEvidence_mcmc_hw_rmatrix", (DL_FUNC) &_graphicalEvidence_mcmc_hw_rmatrix, 9},
    {"_graphicalEvidence_mcmc_last_col", (DL_FUNC) &_graphicalEvidence_mcmc_last_col, 11},
    {"_graphicalEvidence_set_cores", (DL_FUNC) &_graphicalEvidence_set_cores, 1},
    {"_graphicalEvidence_set_seed", (DL_FUNC) &_graphicalEvidence_set_seed, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_graphicalEvidence(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
