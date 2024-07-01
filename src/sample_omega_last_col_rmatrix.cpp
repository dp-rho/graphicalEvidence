#include "graphical_evidence.h"


/*
 * Sample Omega using Hao Wang decomposition Wang decomposition MCMC sampling.
 * Updates current reduced omega using last col restricted sampler, this
 * function operates BGL, GHS, and Wishart priors
 */

void sample_omega_last_col_rmatrix(
  const unsigned int p_reduced,
  const double shape_param,
  const double* scale_params,
  const double omega_pp,
  const double lambda,
  const double dof,
  const int prior,
  arma::vec& beta,
  arma::mat& omega_reduced,
  arma::mat& inv_c,
  arma::mat& inv_omega_11,
  arma::mat& sigma,
  arma::mat& tau,
  arma::mat& nu,
  arma::umat const& ind_noi_mat,
  arma::mat const& s_mat,
  arma::mat const& gibbs_mat,
  arma::mat const& last_col_outer
) {

  /* Tilda parameterization sampling  */
  omega_reduced -= ((1 / omega_pp) * last_col_outer);

  /* Allow selection of all elements besides the ith element  */
  arma::uvec ind_noi;

  /* Needed for GHS case  */
  double lambda_sq = 0;

  /* LAPACK/BLAS args */
  int lapack_dim = p_reduced - 1;

  /* Use existing global memory to avoid constant reallocaiton  */
  arma::vec flex_mem = arma::vec(g_vec1, p_reduced - 1, false, true);
  arma::vec solve_for = arma::vec(g_vec2, p_reduced - 1, false, true);

  /* Iterate through 1 to p_reduced for restricted sampler  */
  for (arma::uword i = 0; i < p_reduced; i++) {

    /* Use existing global memory to avoid constant reallocaiton  */
    ind_noi = ind_noi_mat.unsafe_col(i);

    /* Get sampled gamma value  */
    double gamma_sample = g_rgamma.GetSample(shape_param, scale_params[i]);
    // double gamma_sample = extract_rgamma();

    /* Update inv_omega_11  */
    g_last_col_t1.TimerStart();
    efficient_inv_omega_11_calc(inv_omega_11, ind_noi, sigma, p_reduced, i);
    // inv_omega_11 = arma::inv_sympd(omega_reduced.submat(ind_noi, ind_noi));
    g_last_col_t1.TimerEnd();

    /* Dependent on prior, intialize solve_for in g_vec2*/
    
    /* Wishart case */
    if (prior == WISHART) {

      /* Not implemented  */
    }
    
    /* BGL case */
    else if (prior == BGL) {

      /* Set inv_c */
      inv_c = inv_omega_11 * (lambda + s_mat.at(i ,i));

      /* Calculate lambda^2 one time  */
      lambda_sq = lambda * lambda;

      /* Set solve for  */

      /* Update diagonal of inv_c by adding 1 / tau_12  */
      for (unsigned int j = 0; j < p_reduced - 1; j++) {

        /* Update diagonal of inv_c by adding 1 / tau_12  */
        inv_c.at(j, j) += (1 / tau.at(ind_noi[j], i));

        /* Update b where inv_c x = b */
        solve_for[j] = (
          s_mat.at(ind_noi[j], i) + (
            gibbs_mat.at(ind_noi[j], i) +
            (last_col_outer.at(ind_noi[j], i) / omega_pp)
          ) / tau.at(ind_noi[j], i)
        );
      }
    }

    /* GHS case */
    else if (prior == GHS) {

      /* Set inv_c */
      inv_c = inv_omega_11 * ((1 / lambda) + s_mat.at(i, i));

      /* Calculate lambda^2 one time  */
      lambda_sq = lambda * lambda;

      /* Update inv_c and solve_for */
      for (unsigned int j = 0; j < p_reduced - 1; j++) {

        /* Update diagonal of inv_c by adding 1 / tau_12 * lambda^2 */
        inv_c.at(j, j) += ((1 / tau.at(ind_noi[j], i)) * lambda_sq);

        /* Update b where inv_c x = b */
        solve_for[j] = (
          s_mat.at(ind_noi[j], i) + (
            gibbs_mat.at(ind_noi[j], i) + 
            (last_col_outer.at(ind_noi[j], i) / omega_pp)
          ) / (tau.at(ind_noi[j], i) * lambda_sq)
        );
      }

      // arma::cout << "compiled solve for (mu_i):\n " << solve_for << arma::endl;
    }
    
    /* -mu_i = solve(inv_c, solve_for), store chol(inv_c) in the pointer of inv_c */
    g_last_col_t2.TimerStart();
    LAPACK_dposv(
      &uplo, &lapack_dim, &nrhs, inv_c.memptr(), &lapack_dim, solve_for.memptr(), 
      &lapack_dim, &info_int
    );
    g_last_col_t2.TimerEnd();
    /* g_last_col_t2.TimerStart();
    arma::vec mu_i = arma::solve(inv_c, solve_for);
    g_last_col_t2.TimerEnd(); */

    /* Generate random normals needed to solve for beta */
    flex_mem.randn();
    // extract_rnorm(flex_mem.memptr(), lapack_dim);

    /* Solve chol(inv_c) x = randn(), store result in flex_mem  */
    cblas_dtrsm(
      CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, lapack_dim, 
      nrhs, one, inv_c.memptr(), lapack_dim, flex_mem.memptr(), lapack_dim
    );
    /* g_last_col_t5.TimerStart();
    inv_c = arma::chol(inv_c);
    arma::vec temp = arma::solve(inv_c, flex_mem);
    g_last_col_t5.TimerEnd(); */

    /* Update beta  */
    beta = -solve_for + flex_mem;
    // beta = temp - mu_i;

    /* Update ith col and row of omega and calculate                          */
    /* omega_22 = gamma_sample + (beta.t() * inv_omega_11 * beta) in flex_mem */
    double omega_22 = gamma_sample;
    for (unsigned int j = 0; j < (p_reduced - 1); j++) {

      /* Update the col and row indices excluding the diagonal  */
      omega_reduced.at(ind_noi[j], i) = beta[j];
      omega_reduced.at(i, ind_noi[j]) = beta[j];

      /* Store beta.t() * inv_omega_11 in g_vec1  */
      flex_mem[j] = 0;
      for (unsigned int k = 0; k < (p_reduced - 1); k++) {

        /* First beta.t() * inv_omega_11[, k] */
        flex_mem[j] += (beta[k] * inv_omega_11.at(k, j));
      }

      /* Accumulate beta_omega[j] * beta[j] */
      omega_22 += (beta[j] * flex_mem[j]);
    }
    omega_reduced.at(i, i) = omega_22;

    /* Conditional on prior, update sampling memory */
    
    /* Wishart case */
    if (prior == WISHART) {
      /* TODO: Wishart*/
    }

    /* BGL case */
    else if (prior == BGL) {
      
      /* Calculate a_gig_tau and resuse variable to sample tau_12 */
      for (unsigned int j = 0; j < (p_reduced - 1); j++) {

        double gig_tau = pow(
          beta[j] + gibbs_mat.at(i, ind_noi[j]) + 
          (last_col_outer.at(i, ind_noi[j]) / omega_pp),
          2
        );
        gig_tau = 1 / gigrnd(-1.0 / 2.0, gig_tau, lambda_sq);
        tau.at(ind_noi[j], i) = gig_tau;
        tau.at(i, ind_noi[j]) = gig_tau;
      }
    }

    /* GHS case */
    else if (prior == GHS) {
      
      g_last_col_t6.TimerStart();

      /* Sample tau_12 and nu_12 */
      for (unsigned int j = 0; j < (p_reduced - 1); j++) {

        const double cur_rate = (
          pow(beta[j], 2) / (2 * lambda_sq) + (1 / nu.at(ind_noi[j], i))
        );

        const double cur_tau = 1 / g_rgamma.GetSample(1, 1 / cur_rate);
        // const double cur_tau = 1 / extract_rgamma();
        const double cur_nu = 1 / g_rgamma.GetSample(1, 1 / (1 + (1 / cur_tau)));
        // const double cur_nu = 1 / extract_rgamma();

        tau.at(ind_noi[j], i) = cur_tau;
        tau.at(i, ind_noi[j]) = cur_tau;
        nu.at(ind_noi[j], i) = cur_nu;
        nu.at(i, ind_noi[j]) = cur_nu;
      }

      g_last_col_t6.TimerEnd();
    }

    /* After omega_reduced is updated, update sigma where beta.t() %*% inv_omega_11 */
    /* is still stored in global memory g_vec1 from our previous update of omega    */

    update_sigma_inplace(
      sigma, inv_omega_11, g_vec1, ind_noi, gamma_sample, p_reduced, i
    );
    g_last_col_t1.TimerEnd();

  }

  /* Update omega reduced with last col outer product */
  omega_reduced += ((1 / omega_pp) * last_col_outer);

}