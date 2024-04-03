#include "graphical_evidence.h"

/*
 * Custom memory efficient implementation to update solve_for  
 */

void in_place_update_solve_for(
  arma::rowvec const& rvec,
  arma::mat const& smat,
  arma::uvec const& row_indices,
  arma::uvec const& col_indices
) {

  for (unsigned int i = 0; i < col_indices.n_elem; i++) {
    double sum = 0.0;
    for (unsigned int j = 0; j < row_indices.n_elem; j++) {
      sum += rvec[j] * smat.at(row_indices[j], col_indices[i]);
    }
    g_vec2[i] += sum;
  }
}


/*
 * Custom memory efficient implementation to solve for mu_reduced
 */

void in_place_solve_mu_reduced(
  arma::mat const& inv_c,
  arma::uvec const& one_indices
) {

  int dim = one_indices.n_elem;
  int assign_index = 0;
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      g_mat1[assign_index++] = inv_c.at(one_indices[j], one_indices[i]);
    }
  }

  /* reuse assign index as NRHS */
  assign_index = 1;
  LAPACK_dgesv(&dim, &assign_index, g_mat1, &dim, g_ipv, g_vec2, &dim, &info_int);
}


/*
 * Efficient manual memory management implementation of solving for mu reduced 
 */

void solve_mu_reduced_hw_in_place(
  const unsigned int i,
  arma::uvec const& reduced_one_ind,
  arma::uvec const& reduced_zero_ind,
  arma::mat const& inv_c,
  arma::rowvec const& vec_acc_21
) {

  /* If there are any zero values in currently considered column  */
  if (reduced_zero_ind.n_elem) {

    /* Multiply in place with relevant inv_C submat */
    in_place_update_solve_for(
      vec_acc_21, inv_c, reduced_zero_ind, reduced_one_ind
    );
  }

  /* Solve for mu i reduced */
  in_place_solve_mu_reduced(
    inv_c, reduced_one_ind
  );
}

/*
 * Calculate the mean mu of the multivariate density for eq9 in Hao
 * Wang decomposition MCMC sampling.
 * Native Armadillo implementation
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