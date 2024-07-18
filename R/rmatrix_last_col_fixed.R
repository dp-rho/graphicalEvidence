
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
  
  # Initialize time to calculate restricted Hao Wang sampler
  start_time_mcmc_last_col <- proc.time()
  
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
  
  # ### R code time profiling ###
  calc_time <- proc.time() - start_time_mcmc_last_col
  g_time_env$mcmc_last_col_calc_time <- (
    g_time_env$mcmc_last_col_calc_time + calc_time
  )
  
  cat("last col calc time: \n")
  print(calc_time)

  return(list(
    MC_avg_eq_11=ans_last_col[[1]],
    start_point_first_gibbs=ans_last_col[[2]],
    post_mean_omega_22=ans_last_col[[3]]
  ))
  
}