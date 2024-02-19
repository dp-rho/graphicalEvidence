#include "graphical_evidence.h"

/* Global intializations for compiled timers  */
FunctionTimer g_eq_9_timer("mc_avg_eq_9");


/*
 * Top level function to print currently considered timers
 */
// [[Rcpp::export]]
void print_times(const int nruns) {
  g_eq_9_timer.getTotalDuration(nruns);

}

/*
 * Top level function to reset considered timers
 */
// [[Rcpp::export]]
void reset_times() {
  g_eq_9_timer.resetDuration();

}
