#include "graphical_evidence.h"


/*
 * Sample Omega using Hao Wang decomposition Wang decomposition MCMC sampling.
 * Updates omega accumulator, inv_c accumulator, and mean_vec accumulator as 
 * determined by iterations exceeding burnin period
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
  arma::mat const& s_mat
) {

  /* Number of cols to be iterated through  */
  arma::uword const p = s_mat.n_rows;

  for (arma::uword i = 0; i < p; i++) {

    /* Generate random gamma sample based */
    const double gamma_sample = g_rgamma.GetSample(
      (double) alpha + (n / 2) + 1, 2 / (s_mat.at(i, i) + scale_mat.at(i, i))
    );

    /* Time profiling */
    //g_inv_omega_11_hw.TimerStart();

    /* Calculate inverse of omega excluding row and col i */
    inv_omega_11 = arma::inv(
      omega.submat(ind_noi_mat.col(i), ind_noi_mat.col(i))
    );

    /* Time profiling */
    //g_inv_omega_11_hw.TimerEnd();

    /* Update ith row and col of omega by calculating beta, if i == (p - 1),  */
    /* and burnin is completed, save results in accumulator variables         */

    /* Fill in any zero indices (may be empty)  */
    arma::rowvec vec_acc_21 = -arma::vec(gibbs_mat.submat(
      ind_noi_mat.col(i), arma::uvec({i})
    )).elem(find_which_zeros[i]).t();

    /* Set vec to relevant elements of beta */
    beta.elem(find_which_zeros[i]) = vec_acc_21.t();

    /* Calculate mean mu and assign to one indices if they exist */
    const arma::uword reduced_dim = find_which_ones[i].n_elem;
    if (reduced_dim) {

      /* Time profiling */
      //g_mu_reduced1_hw.TimerStart();

      /* Case where some ones are found in current col  */
      inv_c = inv_omega_11 * (scale_mat.at(i, i) + s_mat.at(i, i));

      /* Manual memory management, initialize solve for vector to find mu reduced */
      for (unsigned int j = 0; j < reduced_dim; j++) {
        int row_index = ind_noi_mat.at(find_which_ones[i][j], i);
        g_vec2[j] = s_mat.at(row_index, i) + scale_mat.at(row_index, i);
      }

      solve_mu_reduced_hw_in_place(
        i, find_which_ones[i], find_which_zeros[i], inv_c, vec_acc_21
      );
      
      /* Time profiling */
      //g_mu_reduced1_hw.TimerEnd();

      /* Time profiling */
      //g_mu_reduced2_hw.TimerStart();

      arma::vec mu_reduced = -arma::vec(g_vec2, reduced_dim, false);

      /* Store results if the column considered (i) is the last (p) */
      if (((iter - burnin) >= 0) && (i == (p - 1))) {
        inv_c_required_store.slice((iter - burnin)) = inv_c.submat(
          find_which_ones[i], find_which_ones[i]
        );
        mean_vec_store.col((iter - burnin)) = mu_reduced;
      }

      /* Time profiling */
      //g_mu_reduced2_hw.TimerEnd();

      /* Time profiling */
      //g_mu_reduced3_hw.TimerStart();
      
      mu_reduced += arma::solve(
        arma::chol(inv_c.submat(find_which_ones[i], find_which_ones[i])), 
        arma::randn<arma::vec>(reduced_dim)
      );

      beta.elem(find_which_ones[i]) = mu_reduced;

      /* Time profiling */
      //g_mu_reduced3_hw.TimerEnd();
    }

    /* Time profiling */
    //g_update_omega_hw.TimerStart();

    /* Update ith col and row of omega */
    omega.submat(arma::uvec({i}), ind_noi_mat.col(i)) = beta.t();
    omega.submat(ind_noi_mat.col(i), arma::uvec({i})) = beta;
    omega.at(i, i) = gamma_sample + arma::dot(beta.t() * inv_omega_11, beta);

    /* Time profiling */
    //g_update_omega_hw.TimerEnd();
  }

  /* If iteration is past burnin period, accumulate Omega */
  if ((iter - burnin) >= 0) {
    omega_save += omega;
  }
  
}