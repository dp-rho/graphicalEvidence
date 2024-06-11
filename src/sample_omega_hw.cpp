#include "graphical_evidence.h"


/*
 * Sample Omega using Hao Wang decomposition Wang decomposition MCMC sampling.
 * Updates omega accumulator, inv_c accumulator, and mean_vec accumulator as 
 * determined by iterations exceeding burnin period, G-Wishart prior
 */

void sample_omega_hw(
  const int iter,
  const int burnin,
  const int n,
  const int alpha,
  arma::vec& beta,
  arma::mat& omega,
  arma::mat& inv_omega_11,
  arma::mat& inv_c,
  arma::mat& omega_save,
  arma::mat& mean_vec_store,
  arma::cube& inv_c_required_store,
  arma::mat const& gibbs_mat,
  arma::mat const& g_mat_adj,
  arma::umat const& ind_noi_mat,
  std::vector<arma::uvec> const& find_which_ones,
  std::vector<arma::uvec> const& find_which_zeros,
  arma::mat const& scale_mat,
  arma::mat const& s_mat,
  arma::mat& sigma
) {

  /* Number of cols to be iterated through  */
  arma::uword const p = s_mat.n_rows;

  /* Allow selection of all elements besides the ith element  */
  arma::uvec ind_noi;

  /* Gamma generation constants */
  /* TODO: Move to top scope */
  double scale_params[p];
  const double shape_param = alpha + ((double)n / 2) + 1;
  for (unsigned int i = 0; i < p; i++) {
    scale_params[i] = 2 / (s_mat.at(i, i) + scale_mat.at(i, i));
  }

  for (arma::uword i = 0; i < p; i++) {

    /* Use existing global memory to avoid constant reallocaiton  */
    ind_noi = ind_noi_mat.unsafe_col(i);

    /* Number of ones in current col */
    const arma::uword reduced_dim = find_which_ones[i].n_elem;

    /* Generate random gamma sample based */
    const double gamma_sample = g_rgamma.GetSample(
      shape_param, scale_params[i]
    );

    /* Time profiling */
    g_inv_omega_11_hw.TimerStart();

    /* Inverse of omega excluding row i and col i can be solved in O(n^2) */
    efficient_inv_omega_11_calc(inv_omega_11, ind_noi, sigma, p, i);

    /* Time profiling */
    g_inv_omega_11_hw.TimerEnd();

    /* Update ith row and col of omega by calculating beta, if i == (p - 1),  */
    /* and burnin is completed, save results in accumulator variables         */

    /* Fill in any zero indices (may be empty) to g_vec1  */
    for (unsigned int j = 0; j < find_which_zeros[i].n_elem; j++) {

      /* Current index of zero  */
      const unsigned int which_zero = find_which_zeros[i][j];

      /* Store vec_acc_21 in g_vec1 */
      const double extracted_gibbs = -gibbs_mat.at(ind_noi[which_zero], i);
      g_vec1[j] = extracted_gibbs;

      /* Set vec to relevant elements of beta */
      beta[which_zero] = extracted_gibbs;
    }

    /* Calculate mean mu and assign to one indices if they exist */
    if (reduced_dim) {

      /* Time profiling */
      g_mu_reduced1_hw.TimerStart();

      /* Case where some ones are found in current col  */
      inv_c = inv_omega_11 * (scale_mat.at(i, i) + s_mat.at(i, i));

      /* Manual memory management, initialize solve for vector to find mu reduced */
      for (unsigned int j = 0; j < reduced_dim; j++) {
        int row_index = ind_noi_mat.at(find_which_ones[i][j], i);
        g_vec2[j] = s_mat.at(row_index, i) + scale_mat.at(row_index, i);
      }

      solve_mu_reduced_hw_in_place(
        i, find_which_ones[i], find_which_zeros[i], inv_c
      );

      /* Time profiling */
      g_mu_reduced1_hw.TimerEnd();

      /* Time profiling */
      g_mu_reduced2_hw.TimerStart();

      /* Store results if the column considered (i) is the last (p) */
      if (((iter - burnin) >= 0) && (i == (p - 1))) {

        inv_c_required_store.slice((iter - burnin)) = inv_c.submat(
          find_which_ones[i], find_which_ones[i]
        );

        double* cur_col = mean_vec_store.colptr(iter - burnin);
        memcpy(cur_col, g_vec2, reduced_dim * sizeof(double));
      }

      /* Time profiling */
      g_mu_reduced2_hw.TimerEnd();

      /* Time profiling */
      g_mu_reduced3_hw.TimerStart();

      /* Generate random normals in g_vec1  */
      for (unsigned int j = 0; j < reduced_dim; j++) {
        g_vec1[j] = arma::randn();
      }

      /* Solve chol(inv_c) x = randn(), store result in g_vec1  */
      cblas_dtrsm(
        CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, reduced_dim, nrhs, one,
        g_mat1, reduced_dim, g_vec1, reduced_dim
      );

      /* Update one indices of beta with mu_i + solve(chol(inv_c_ones, randn()))  */
      for (unsigned int j = 0; j < reduced_dim; j++) {
        beta[find_which_ones[i][j]] = g_vec2[j] + g_vec1[j];
      }

      /* Time profiling */
      g_mu_reduced3_hw.TimerEnd();
    }

    /* Time profiling */
    g_update_omega_hw.TimerStart();

    /* Update ith col and row of omega and */
    /* calculate omega_22 = gamma_sample + (beta.t() * inv_omega_11 * beta) in g_vec1 */
    double omega_22 = gamma_sample;
    for (unsigned int j = 0; j < (p - 1); j++) {

      /* Update the col and row indices excluding the diagonal  */
      omega.at(ind_noi[j], i) = beta[j];
      omega.at(i, ind_noi[j]) = beta[j];

      /* Store beta.t() * inv_omega_11 in g_vec1  */
      g_vec1[j] = 0;
      for (unsigned int k = 0; k < (p - 1); k++) {

        /* First beta.t() * inv_omega_11[, k] */
        g_vec1[j] += (beta[k] * inv_omega_11.at(k, j));
      }

      /* Accumulate beta_omega[j] * beta[j] */
      omega_22 += (beta[j] * g_vec1[j]);
    }
    omega.at(i, i) = omega_22;

    update_sigma_inplace(
      sigma, inv_omega_11, g_vec1, ind_noi, gamma_sample, p, i
    );

    /* Time profiling */
    g_update_omega_hw.TimerEnd();
  }

  /* If iteration is past burnin period, accumulate Omega */
  if ((iter - burnin) >= 0) {
    omega_save += omega;
  }
  
}