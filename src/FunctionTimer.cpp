#include "graphical_evidence.h"

/* Global intializations for compiled timers  */
FunctionTimer g_eq_9_timer("mc_avg_eq_9");
FunctionTimer g_mcmc_hw_timer("mcmc_hw");
FunctionTimer g_eq_11_timer("mc_avg_eq_11");
FunctionTimer g_mcmc_last_col_timer("mcmc_last_col");
FunctionTimer g_mu_reduced1_hw("mu_reduced1_hw");
FunctionTimer g_mu_reduced2_hw("mu_reduced2_hw");
FunctionTimer g_mu_reduced3_hw("mu_reduced3_hw");
FunctionTimer g_update_omega_hw1("update_omega_hw1");
FunctionTimer g_update_omega_hw2("update_omega_hw2");
FunctionTimer g_sample_omega_hw("sample_omega_hw");
FunctionTimer g_inv_omega_11_hw("inv_omega_11_hw");
FunctionTimer g_inv_c_hw("g_inv_c_hw");
FunctionTimer g_last_col_t1("g_last_col_t1");
FunctionTimer g_last_col_t2("g_last_col_t2");
FunctionTimer g_last_col_t3("g_last_col_t3");
FunctionTimer g_last_col_t4("g_last_col_t4");
FunctionTimer g_last_col_t5("g_last_col_t5");
FunctionTimer g_last_col_t6("g_last_col_t6");


/*
 * Top level function to print currently considered timers
 */
// [[Rcpp::export]]
void print_times(const int nruns) {
  g_eq_9_timer.getTotalDuration(nruns);
  g_mcmc_hw_timer.getTotalDuration(nruns);
  g_eq_11_timer.getTotalDuration(nruns);
  g_mcmc_last_col_timer.getTotalDuration(nruns);
  g_mu_reduced1_hw.getTotalDuration(nruns);
  g_mu_reduced2_hw.getTotalDuration(nruns);
  g_mu_reduced3_hw.getTotalDuration(nruns);
  g_update_omega_hw1.getTotalDuration(nruns);
  g_update_omega_hw2.getTotalDuration(nruns);
  g_sample_omega_hw.getTotalDuration(nruns);
  g_inv_omega_11_hw.getTotalDuration(nruns);
  g_inv_c_hw.getTotalDuration(nruns);
  g_last_col_t1.getTotalDuration(nruns);
  g_last_col_t2.getTotalDuration(nruns);
  g_last_col_t3.getTotalDuration(nruns);
  g_last_col_t4.getTotalDuration(nruns);
  g_last_col_t5.getTotalDuration(nruns);
  g_last_col_t6.getTotalDuration(nruns);
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
  g_mu_reduced1_hw.resetDuration();
  g_mu_reduced2_hw.resetDuration();
  g_mu_reduced3_hw.resetDuration();
  g_update_omega_hw1.resetDuration();
  g_update_omega_hw2.resetDuration();
  g_sample_omega_hw.resetDuration();
  g_inv_omega_11_hw.resetDuration();
  g_inv_c_hw.resetDuration();
  g_last_col_t1.resetDuration();
  g_last_col_t2.resetDuration();
  g_last_col_t3.resetDuration();
  g_last_col_t4.resetDuration();
  g_last_col_t5.resetDuration();
  g_last_col_t6.resetDuration();
}
