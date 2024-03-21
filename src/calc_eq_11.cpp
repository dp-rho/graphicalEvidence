#include "graphical_evidence.h"

/*
 * Calculates equation 11 in paper using sampled results for omega_reduced
 * and the elements of omega at column and row p throughout the MCMC sampling.
 * This can be parallelized with OpenMP as the calculation is not inherently
 * sequential.
 */
// [[Rcpp::export]]
double calc_eq_11(
  const double omega_22_mean,
  const double s_22,
  const double scale_mat_22,
  const double alpha,
  const unsigned int n,
  const unsigned int nmc,
  NumericVector gamma_vec
) {

  // ARG
  // arma::vec const& gamma_subtractors

  /* Create armadillo constructor object  */
  arma::vec gamma_subtractors(gamma_vec.begin(), nmc);

  /* Time profiling */
  g_eq_11_timer.TimerStart();
  
  /* Shape and scale parameter for gamma density  */
  const double scale_param = 2 / (s_22 + scale_mat_22);
  const double shape_param = alpha + ((double) n / 2) + 1;

  /* Iterate through accumulated gamma subtractors and calculate gamma density  */
  double gamma_acc = 0.0;
  for (arma::uword i = 0; i < nmc; i++) {

    /* Gamma value to check is mean pth row/col of omega minus saved subtractors  */
    double temp_gamma = omega_22_mean - gamma_subtractors[i];

    /* Check if temp_gamma is not positive, if so, the technique is not applicable  */

    /* Accumulate the gamma density */
    gamma_acc += gamma_density(temp_gamma, shape_param, scale_param);
  }
  gamma_acc /= nmc;

  /* Time profiling */
  g_eq_11_timer.TimerEnd();

  return log(gamma_acc);
}