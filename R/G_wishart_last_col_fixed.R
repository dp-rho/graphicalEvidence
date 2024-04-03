
G_wishart_last_col_fixed <- function(
  S,
  n,
  burnin,
  nmc,
  alpha,
  fixed_last_col,
  scale_matrix,
  G_mat_adj,
  matrix_accumulator_gibbs,
  post_mean_omega
) {
  
  # Initialize time to calculate restricted Hao Wang sampler
  start_time_mcmc_last_col <- proc.time()
  
  p <- nrow(S)
  
  ##################################
  # RcppArmadillo implementation ###
  ##################################
  
  # RcppArmadillo implementation
  ans_last_col <- mcmc_last_col(
    n, burnin, nmc, alpha, p, fixed_last_col, S, scale_matrix, G_mat_adj,
    matrix_accumulator_gibbs, post_mean_omega
  )
  
  ##################################
  
  ### R code time profiling ###
  calc_time <- proc.time() - start_time_mcmc_last_col
  g_time_env$mcmc_last_col_calc_time <- (
    g_time_env$mcmc_last_col_calc_time + calc_time
  )
  ######################
  
  return(list(
    MC_average_Equation_11=ans_last_col[[1]],
    start_point_first_gibbs=ans_last_col[[2]],
    post_mean_omega_22=ans_last_col[[3]]
  ))
}