#include "graphical_evidence.h"


/*
 * Sample Omega using Hao Wang decomposition Wang decomposition MCMC sampling.
 * Updates current reduced omega using last col restricted sampler
 */

void sample_omega_last_col(
  const unsigned int p_reduced,
  const double shape_param,
  const double* scale_params,
  const double omega_pp,
  arma::vec& beta,
  arma::mat& omega_reduced,
  arma::mat& inv_c,
  arma::vec& ind_noi_data1,
  arma::vec& ind_noi_data2,
  std::vector<arma::uvec> const& find_which_ones,
  std::vector<arma::uvec> const& find_which_zeros,
  arma::umat const& ind_noi_mat,
  arma::mat const& s_mat,
  arma::mat const& scale_mat,
  arma::mat const& gibbs_mat,
  arma::mat const& last_col_outer
) {

  omega_reduced += ((1 / omega_pp) * last_col_outer);

  /* Data vec for reduced (which ones) in sampler */
  arma::vec reduced_vec;
  arma::mat inv_c_required;

  /* Iterate through 1 to p_reduced for restricted sampler  */
  for (arma::uword i = 0; i < p_reduced; i++) {

    /* Extract no i indices for ith col */
    ind_noi_data1 = -(1 / omega_pp) * last_col_outer.submat(
      ind_noi_mat.col(i), arma::uvec({i})
    );
    ind_noi_data2 = gibbs_mat.submat(ind_noi_mat.col(i), arma::uvec({i}));
    double gamma_sample = g_rgamma.GetSample(shape_param, scale_params[i]);
    // double gamma_sample = extract_rgamma();

    /* Update inv_c */
    inv_c = arma::inv(
      omega_reduced.submat(ind_noi_mat.col(i), ind_noi_mat.col(i))
    ) * s_mat.at(i, i) * scale_mat.at(i, i);

    const arma::uword reduced_dim = find_which_ones[i].n_elem;
    if (reduced_dim) {

      inv_c_required = inv_c.submat(find_which_ones[i], find_which_ones[i]);

      if (find_which_zeros[i].n_elem) {

        /* Update zero indices in beta  */
        beta.elem(find_which_zeros[i]) = -1 * (
          ind_noi_data1.elem(find_which_zeros[i]) +
          ind_noi_data2.elem(find_which_zeros[i])
        );

        arma::mat inv_c_not_required = inv_c.submat(
          find_which_zeros[i], find_which_ones[i]
        );

        reduced_vec = -solve(
          inv_c_required,
          (
            arma::vec(scale_mat.submat(ind_noi_mat.col(i), arma::uvec({i}))).elem(
              find_which_ones[i]
            ) +
            arma::vec(s_mat.submat(ind_noi_mat.col(i), arma::uvec({i}))).elem(
              find_which_ones[i]
            ) +
            (ind_noi_data2.elem(find_which_zeros[i]).t() * inv_c_not_required).t() +
            (ind_noi_data1.elem(find_which_zeros[i]).t() * inv_c_not_required).t()
          )
        );
      }
      else {
        //arma::vec rnorm_vec = arma::zeros(reduced_dim);
        //extract_rnorm(rnorm_vec.memptr(), reduced_dim);
        reduced_vec = -solve(
          inv_c_required,
          arma::randn<arma::vec>(reduced_dim)
        );
      }
      
      /* Update beta at indices associated with 1's */
      //arma::vec rnorm_vec = arma::zeros(reduced_dim);
      //extract_rnorm(rnorm_vec.memptr(), reduced_dim);
      reduced_vec += arma::solve(
        arma::chol(inv_c_required), arma::randn<arma::vec>(reduced_dim)
      );
      beta.elem(find_which_ones[i]) = reduced_vec;

    }
    else {
      beta = (
        (-1 * gibbs_mat.submat(ind_noi_mat.col(i), arma::uvec({i}))) +
        (-1 / omega_pp * last_col_outer.submat(
          ind_noi_mat.col(i), arma::uvec({i})
        ))
      );
    }

    /* Update row and col i using sampled values  */
    omega_reduced.submat(ind_noi_mat.col(i), arma::uvec({i})) = beta;
    omega_reduced.submat(arma::uvec({i}), ind_noi_mat.col(i)) = beta.t();
    omega_reduced.at(i, i) = gamma_sample + arma::dot(beta.t(),
      inv_c / (s_mat.at(i, i) * scale_mat.at(i, i)) * beta
    );
  }

  /* Update omega reduced with last col outer product */ 
  omega_reduced += ((1 / omega_pp) * last_col_outer);
}