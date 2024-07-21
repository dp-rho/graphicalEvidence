
# Last col restricted sampler for BGL
rmatrix_last_col_fixed <- function(
  S,
  n,
  burnin,
  nmc,
  fixed_last_col,
  prior,
  dof = 0,
  matrix_acc_gibbs = 0,
  post_mean_omega = 0,
  post_mean_tau = 0,
  lambda = 0
) {
  
  # Get arguments for C++ call
  p <- nrow(S)
  
  coded_prior <- switch(
    prior,
    'Wishart' = 0,
    'BGL' = 1,
    'GHS'= 2
  )
  
  ##################################
  # RcppArmadillo implementation  ##
  ##################################

  # RcppArmadillo implementation
  ans_last_col <- mcmc_last_col_rmatrix(
    n, burnin, nmc, p, dof, lambda, coded_prior, fixed_last_col, S,
    post_mean_tau, matrix_acc_gibbs, post_mean_omega
  )

  ##################################

  return(list(
    MC_avg_eq_11=ans_last_col[[1]],
    start_point_first_gibbs=ans_last_col[[2]],
    post_mean_omega_22=ans_last_col[[3]]
  ))
  
}