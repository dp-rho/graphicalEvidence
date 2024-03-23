#include "graphical_evidence.h"


/*
 * Use Hao Wang decomposition to run an MCMC sampler nmc + burnin times
 * and accumulate omega, mean_vec, and inv_c_required to calculate
 * mc average equation 9 after mcmc sampling is complete.
 * Currently communicates with R through Rcpp interface, will
 * eventually communicate only with other compiled code
 */
 // [[Rcpp::export]]
List mcmc_last_col(
  const unsigned int n,
  const unsigned int burnin,
  const unsigned int nmc,
  const unsigned int alpha,
  const unsigned int p,
  NumericVector fixed_last_col_nvec,
  NumericVector s_mat_nvec,
  NumericVector scale_mat_nvec,
  NumericVector g_mat_adj_nvec,
  NumericVector gibbs_mat_nvec,
  NumericVector post_mean_omega_nvec
) {

  /* Deep copy Rcpp objects to Armadillo constructs */
  arma::vec fixed_last_col(fixed_last_col_nvec.begin(), p - 1);
  arma::mat s_mat(s_mat_nvec.begin(), p, p);
  arma::mat scale_mat(scale_mat_nvec.begin(), p, p);
  arma::mat omega(post_mean_omega_nvec.begin(), p, p);
  arma::mat gibbs_mat(gibbs_mat_nvec.begin(), p, p);
  arma::mat g_mat_adj(g_mat_adj_nvec.begin(), p, p);

  /* Time profiling */
  g_mcmc_last_col_timer.TimerStart();

  /* Set omega_22 gamma parameters  */
  double scale_params[p];
  const double shape_param = alpha + ((double) n / 2) + 1;
  for (unsigned int i = 0; i < p; i++) {
    scale_params[i] = 2 / (s_mat.at(i, i) + scale_mat.at(i, i));
  }

  /* Calculate outer product of last col  */
  arma::mat last_col_outer = fixed_last_col * fixed_last_col.t();

  /* Reduce matrices */
  const int p_reduced = p - 1;
  omega = omega.submat(0, 0, p_reduced - 1, p_reduced - 1);
  gibbs_mat = gibbs_mat.submat(0, 0, p_reduced - 1, p_reduced - 1);
  g_mat_adj = g_mat_adj.submat(0, 0, p_reduced - 1, p_reduced - 1);
  scale_mat = scale_mat.submat(0, 0, p_reduced - 1, p_reduced - 1);
  s_mat = s_mat.submat(0, 0, p_reduced - 1, p_reduced - 1);

  /* Identify indices where ones and zeros exist for adjacency matrix */
  arma::umat ind_noi_mat(p_reduced - 1, p_reduced);
  std::vector<arma::uvec> find_which_ones(p_reduced);
  std::vector<arma::uvec> find_which_zeros(p_reduced);
  initialize_indices(
    g_mat_adj, ind_noi_mat, find_which_ones, find_which_zeros
  );

  /* Initialize accumulators  */
  arma::mat omega_reduced_acc = arma::zeros(p_reduced, p_reduced);
  arma::vec gamma_subtractors = arma::zeros(nmc);
  double omega_22_acc = 0;

  /* Initialize calculation memory  */
  arma::mat inv_c = arma::zeros(p_reduced - 1, p_reduced - 1);
  arma::vec ind_noi_data1 = arma::zeros(p_reduced - 1);
  arma::vec ind_noi_data2 = arma::zeros(p_reduced - 1);
  arma::vec beta = arma::zeros(p_reduced - 1);

  /* Perform main loop of mcmc sampling */
  for (arma::uword i = 0; i < (burnin + nmc); i++) {

    /* Get random gamma value and calculate omega_22  */
    double gamma_sample = g_rgamma.GetSample(shape_param, scale_params[p - 1]);
    // double gamma_sample = extract_rgamma();
    double omega_22 = gamma_sample + arma::dot(
      fixed_last_col.t(), arma::solve(omega, fixed_last_col)
    );

    /* In general case call omega reduced sampler */
    if (p_reduced > 1) {
      sample_omega_last_col(
        p_reduced, shape_param, scale_params, omega_22, beta, omega, inv_c,
        ind_noi_data1, ind_noi_data2, find_which_ones, find_which_zeros,
        ind_noi_mat, s_mat, scale_mat, gibbs_mat, last_col_outer
      );

    }

    /* In special case where the matrices considered are single elements  */
    /* the sampling process is highly simplified                          */
    else {
      gamma_sample = g_rgamma.GetSample(shape_param, scale_params[0]);
      omega.at(0, 0) = gamma_sample + (
        fixed_last_col[0] * fixed_last_col[0] / omega_22
      );
    }

    /* Save results if burnin is complete */
    if (i >= burnin) {
      omega_22_acc += omega_22;
      omega_reduced_acc += omega;
      gamma_subtractors[i - burnin] = arma::dot(
        fixed_last_col.t(), arma::solve(omega, fixed_last_col)
      );
    }
  }

  /* Get posterior means of sampled values */
  omega_22_acc /= nmc;
  omega_reduced_acc /= nmc;

  List z = List::create(
    Rcpp::wrap(gamma_subtractors), Rcpp::wrap(omega_reduced_acc), omega_22_acc
  );

  /* Time profiling */
  g_mcmc_last_col_timer.TimerEnd();

  return z;

}