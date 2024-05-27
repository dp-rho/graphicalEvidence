#include "graphical_evidence.h"


/*
 * Sample Omega using Hao Wang decomposition Wang decomposition MCMC sampling.
 * Updates omega accumulator, inv_c accumulator, and mean_vec accumulator as
 * determined by iterations exceeding burnin period for Wishart, BGL, and
 * GHS prior, with conditional branching determined by the input prior.
 */

void sample_omega_hw_rmatrix(
  const int iter,
  const int burnin,
  const int n,
  const int prior,
  const int dof,
  const int lambda,
  arma::vec& beta,
  arma::mat& omega,
  arma::mat& cur_sigma,
  arma::mat& tau,
  arma::mat& nu,
  arma::mat& inv_omega_11,
  arma::mat& inv_c,
  arma::mat& omega_save,
  arma::mat& tau_save,
  arma::mat& mean_vec_store,
  arma::cube& inv_c_required_store,
  arma::mat const& gibbs_mat,
  arma::umat const& ind_noi_mat,
  arma::mat const& s_mat,
  const double shape,
  const double* scale_vec
) {

  /* Number of cols to be iterated through  */
  arma::uword const p = s_mat.n_rows;

  /* Allow selection of all elements besides the ith element  */
  arma::uvec ind_noi;

  /* Use existing global memory to avoid constant reallocaiton  */
  arma::vec flex_mem = arma::vec(g_vec1, p - 1, false, true);
  arma::vec solve_for = arma::vec(g_vec2, p - 1, false, true);

  /* Iterate p times and update col p and row p at each iteration */
  for (unsigned int i = 0; i < p; i++) {

    /* Use existing global memory to avoid constant reallocaiton  */
    ind_noi = ind_noi_mat.unsafe_col(i);

    /* Update preparation step is dependent on prior  */

    /* Sample from prior specific gamma density */
    const double gamma_param = g_rgamma.GetSample(shape, scale_vec[i]);
    // const double gamma_param = extract_rgamma();

    /* Wishart case */
    if (prior == WISHART) {

      /* Inverse of omega excluding row i and col i is solved directly  */
      inv_omega_11 = arma::inv(omega.submat(ind_noi, ind_noi));
      inv_c = (s_mat.at(i, i) + 1) * inv_omega_11;
      
      /* Solve for vector is just S_21 in Wishart  */
      for (unsigned int j = 0; j < (p - 1); j++) {
        solve_for[j] = s_mat.at(i, ind_noi[j]);
      }
      
    }

    /* BGL case */
    else if (prior == BGL) {

      /* Time profiling */
      g_inv_omega_11_hw.TimerStart();

      /* Inverse of omega excluding row i and col i can be solved in O(n^2) */
      efficient_inv_omega_11_calc(inv_omega_11, ind_noi, cur_sigma, p, i);

      g_inv_omega_11_hw.TimerEnd();

      inv_c = inv_omega_11 * (s_mat.at(i, i) + lambda);
      for (unsigned int j = 0; j < (p - 1); j++) {

        /* Update inv_c diagonal  */
        inv_c.at(j, j) += (1 / tau.at(ind_noi[j], i));

        /* update solve for */
        solve_for[j] = s_mat.at(ind_noi[j], i) + (
          gibbs_mat.at(ind_noi[j], i) / tau.at(ind_noi[j], i)
        );
      }
    }

    /* GHS case */
    else if (prior == GHS) {
      /* TODO: GHS  */
    }

    g_mu_reduced1_hw.TimerStart();
    /* Solve for mu_i, store in flex_mem */
    flex_mem = -arma::solve(inv_c, solve_for);
    g_mu_reduced1_hw.TimerEnd();

    /* Save mu_i and inv_c if this is the final iteration */
    if (((iter - burnin) >= 0) && (i == (p - 1))) {
      inv_c_required_store.slice(iter - burnin) = inv_c;
      mean_vec_store.col(iter - burnin) = flex_mem;
    }

    /* Use existing memory to calculate beta, previous values */
    /* will not be needed again                               */
    g_mu_reduced2_hw.TimerStart();
    inv_c = arma::chol(inv_c);
    solve_for.randn();
    // extract_rnorm(solve_for.memptr(), p - 1);
    beta = flex_mem + arma::solve(inv_c, solve_for);
    g_mu_reduced2_hw.TimerEnd();

    /* Update flex_mem to store inv_omega_11 %*% beta */
    flex_mem = inv_omega_11 * beta;

    /* Update ith row and col of omega  */
    for (unsigned int j = 0; j < (p - 1); j++) {
      omega.at(ind_noi[j], i) = beta[j];
      omega.at(i, ind_noi[j]) = beta[j];
    }
    omega.at(i, i) = gamma_param + arma::dot(beta, flex_mem);
    // arma::cout << "omega at iter: (" << iter << ") ith: (" << i << ")\n" << omega << arma::endl;

    /* Update sampling parameters for next iteration in prior specific method */
    
    /* Wishart case */
    if (prior == WISHART) {
      /* TODO: Wishart */
    }

    /* BGL case */
    else if (prior == BGL) {
      g_mu_reduced3_hw.TimerStart();
      /* Update sigma */
      update_sigma_inplace(
        cur_sigma, inv_omega_11, flex_mem, ind_noi, gamma_param, p, i
      );
      g_mu_reduced3_hw.TimerEnd();
      

      g_update_omega_hw.TimerStart();

      /* Calculate a_gig_tau and resuse variable to sample tau_12 */
      for (unsigned int j = 0; j < (p - 1); j++) {
        //arma::cout << "ith: " << i << " iter " << j << arma::endl;
        double gig_tau = pow(beta[j] + gibbs_mat.at(i, ind_noi[j]), 2);
        //arma::cout << "gig_tau " << gig_tau << arma::endl;
        gig_tau = 1 / gigrnd(-1.0 / 2.0, gig_tau, pow(lambda, 2));
        // arma::cout << "sampled tau: " << gig_tau << arma::endl;
        // gig_tau = 1 / extract_rgig();
        tau.at(ind_noi[j], i) = gig_tau;
        tau.at(i, ind_noi[j]) = gig_tau;
      }

      g_update_omega_hw.TimerEnd();
      
      //arma::cout << "tau at iter: (" << iter << ") ith: (" << i << ")\n" << tau << arma::endl;
    }

    /* GHS case */
    else if (prior == GHS) {
      /* TODO: GHS  */
    }

    // arma::cout << "sigma at iter: (" << iter << ") ith: (" << i << ")\n" << cur_sigma << arma::endl;
    // arma::cout << "tau at iter: (" << iter << ") ith: (" << i << ")\n" << tau << arma::endl;

  }

  /* If iteration is past burnin period, accumulate Omega */
  if ((iter - burnin) >= 0) {
    omega_save += omega;
    tau_save += tau;
  }

}


/* 
 * Computationally efficient calculation of inverse of omega_11:
 * For this implementation, inverse of omega_11 can be calculated as
 * sigma_11 - sigma_12 %*% t(sigma_12) / sigma_22
 */

void efficient_inv_omega_11_calc(
  arma::mat& inv_omega_11,
  arma::uvec const& ind_noi,
  arma::mat const& sigma,
  const unsigned int p,
  const unsigned int ith
) {                

  /* Iterate over cols  */
  for (unsigned int i = 0; i < (p - 1); i++) {

    /* Iterate over rows  */
    for (unsigned int j = 0; j < (p - 1); j++) {

      inv_omega_11.at(j, i) = sigma.at(ind_noi[j], ind_noi[i]) - (
        sigma.at(ind_noi[j], ith) * sigma.at(ind_noi[i], ith) /
        sigma.at(ith, ith)
      );
    }
  }
}


/*
 * Update sigma for either BGL or GHS prior in place using:
 * sigma_11 = inv_omega_11 + (omega_beta %*% t(omega_beta) / gamma_param)
 * sigma_12 = -omega_beta / gamma_param
 * sigma_22 = 1 / gamma_param
 */

void update_sigma_inplace(
  arma::mat& sigma,
  arma::mat const& inv_omega_11,
  arma::vec const& omega_beta,
  arma::uvec const& ind_noi,
  const double gamma_param,
  const unsigned int p,
  const unsigned int ith
) {

  /* Iterate over cols  */
  for (unsigned int i = 0; i < (p - 1); i++) {

    /* Iterate over rows  */
    for (unsigned int j = 0; j < (p - 1); j++) {

      /* Update indices in sigma_11 */
      sigma.at(ind_noi[j], ind_noi[i]) = inv_omega_11.at(j, i) + (
        omega_beta[j] * omega_beta[i] / gamma_param
      );
    }
    
    /* Update indices in sigma_12 */
    sigma.at(ind_noi[i], ith) = -omega_beta[i] / gamma_param;
    sigma.at(ith, ind_noi[i]) = -omega_beta[i] / gamma_param;
  }

  /* Update sigma_22  */
  sigma.at(ith, ith) = 1 / gamma_param;
}