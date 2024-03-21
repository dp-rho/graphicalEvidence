#include "graphical_evidence.h"

/* Global intializations for compiled timers  */
FunctionTimer g_eq_9_timer("mc_avg_eq_9");
FunctionTimer g_mcmc_hw_timer("mcmc_hw");
FunctionTimer g_eq_11_timer("mc_avg_eq_11");
FunctionTimer g_mcmc_last_col_timer("mcmc_last_col");


/*
 * Top level function to print currently considered timers
 */
// [[Rcpp::export]]
void print_times(const int nruns) {
  g_eq_9_timer.getTotalDuration(nruns);
  g_mcmc_hw_timer.getTotalDuration(nruns);
  g_eq_11_timer.getTotalDuration(nruns);
  g_mcmc_last_col_timer.getTotalDuration(nruns);
}


/*
 * Top level function to reset considered timers
 */
// [[Rcpp::export]]
void reset_times() {
  g_eq_9_timer.resetDuration();
  g_mcmc_hw_timer.resetDuration();
  g_eq_11_timer.resetDuration();
  g_mcmc_last_col_timer.resetDuration();
}
