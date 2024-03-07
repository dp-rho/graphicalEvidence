#include "graphical_evidence.h"


/*
 * Calculate the mean mu of the multivariate density for eq9 in Hao
 * Wang decomposition MCMC sampling
 */

arma::vec solve_mu_reduced_hw(
  const unsigned int i,
  arma::uvec const& reduced_one_ind,
  arma::uvec const& reduced_zero_ind,
  arma::mat const& inv_c,
  arma::rowvec& vec_acc_21,
  arma::vec& solve_for
) {

  /* If there are any zero values in currently considered column  */
  if (reduced_zero_ind.n_elem) {

    /* Multiply in place with relevant inv_C submat */
    vec_acc_21 *= inv_c.submat(reduced_zero_ind, reduced_one_ind);

    solve_for += vec_acc_21.t();
  }

  /* Solve for mu i reduced */
  arma::vec mu_reduced = -arma::solve(
    inv_c.submat(reduced_one_ind, reduced_one_ind), solve_for
  );

  return mu_reduced;
}