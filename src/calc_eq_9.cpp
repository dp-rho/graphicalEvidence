#include "graphical_evidence.h"

/*
 * Top level function to calculate MC average of equation 9 in graphical evidence paper
 * This function currently communicates with R and deep copies data to RcppArmadillo
 * constructors, in fully compiled code the deep copy step will not be necessary
 * as this function will be communicating with other compiled code.
 */
// [[Rcpp::export]]
List calc_eq_9(
  NumericVector ones_vec,
  NumericVector one_indices,
  NumericVector post_mean_omega,
  NumericVector inv_sigma_stores,
  NumericVector mean_vecs,
  int p,
  int nmc
 ) {

  /* Deep copy R objects to Armadillo constructors */
  arma::rowvec arma_ones(ones_vec.begin(), p - 1);
  int xdim = arma::sum(arma_ones);
  arma::cube arma_inv_sigma_stores(inv_sigma_stores.begin(), xdim, xdim, nmc);
  arma::mat arma_mean_vecs(mean_vecs.begin(), nmc, xdim);

  /* Time profiling */
  g_eq_9_timer.TimerStart();

  /* If no nonzero entries of adj matrix last col, we can immediately return 0  */
  double mc_avg_eq_9 = 0;
  if (xdim) {
    
    /* Create x data to be evaluated with multivariate normal pdf */
    arma::rowvec x(xdim);

    /* Extract pth row for indices identified as 1's in last col of adj matrix  */
    for (int i = 0; i < xdim; i++) {
      x(i) = post_mean_omega[((one_indices[i] - 1) * p) + (p - 1)];
    }

    /* Calculate log mean of multivariate normal density  */
    mc_avg_eq_9 = log(arma::mean(exp(log_dmvnrm_arma_vec(x, arma_mean_vecs, arma_inv_sigma_stores))));
  }

  List z = List::create(post_mean_omega, mc_avg_eq_9);

  /* Time profiling */
  g_eq_9_timer.TimerEnd();

  return z;
}
