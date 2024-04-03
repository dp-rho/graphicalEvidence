
G_wishart_Hao_wang <- function(
  S,
  n,
  burnin,
  nmc,
  alpha,
  scale_matrix,
  G_mat_adj,
  matrix_accumulator_gibbs,
  start_point_first_gibbs
) {
  
  # Initialize time to calculate Hao Wang sampler
  start_time_mcmc_hw <- proc.time()
  
  # Require matrix type
  S <- as.matrix(S)
  scale_matrix <- as.matrix(scale_matrix)
  G_mat_adj <- as.matrix(G_mat_adj)
  matrix_accumulator_gibbs <- as.matrix(matrix_accumulator_gibbs)
  start_point_first_gibbs <- as.matrix(start_point_first_gibbs)
  
  # Dimension of this iteration of Hao Wang
  p <- nrow(S)
  
  # ##################################
  # # RcppArmadillo implementation ###
  # ##################################
  
  # RcppArmadillo implementation
  ans_hw <- mcmc_hw(
    n, burnin, nmc, alpha, p, S, scale_matrix, G_mat_adj,
    matrix_accumulator_gibbs, start_point_first_gibbs
  )

  ### R code time profiling ###
  calc_time <- proc.time() - start_time_mcmc_hw
  g_time_env$mcmc_hw_calc_time <- (
    g_time_env$mcmc_hw_calc_time + calc_time
  )
  ######################

  return(list(post_mean_omega=ans_hw[[1]], MC_average_Equation_9=ans_hw[[2]]))
}