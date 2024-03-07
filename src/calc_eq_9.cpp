#include "graphical_evidence.h"

/*
 * Top level function to calculate MC average of equation 9 in graphical evidence paper
 */

List calc_eq_9(
  arma::uvec const& find_which_ones,
  arma::cube const& arma_inv_sigma_stores,
  arma::mat const& arma_mean_vecs,
  arma::mat const& arma_mean_omega,
  arma::uword const xdim,
  const unsigned int p,
  const unsigned int nmc
 ) {

  // arma::cout << "mean vecs: " << arma::endl << arma_mean_vecs << arma::endl;
  // arma::cout << "inv sigma stores vecs: " << arma::endl << arma_inv_sigma_stores << arma::endl;

  /* Time profiling */
  g_eq_9_timer.TimerStart();

  /* If no nonzero entries of adj matrix last col, we can immediately return 0  */
  double mc_avg_eq_9 = 0;
  if (xdim) {

    /* Create x data to be evaluated with multivariate normal pdf */
    arma::rowvec x = arma_mean_omega.submat(arma::uvec({p - 1}), find_which_ones);

    /* Calculate log mean of multivariate normal density  */
    mc_avg_eq_9 = log(arma::mean(exp(log_dmvnrm_arma_vec(x, arma_mean_vecs, arma_inv_sigma_stores))));
  }

  List z = List::create(Rcpp::wrap(arma_mean_omega), mc_avg_eq_9);

  /* Time profiling */
  g_eq_9_timer.TimerEnd();

  return z;
}
