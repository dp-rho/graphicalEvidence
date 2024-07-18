
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
  
  # Dimension of this iteration
  p <- nrow(S)
  
  # ##################################
  # # RcppArmadillo implementation ###
  # ##################################
  
  ans_hw <- mcmc_hw(
    n, burnin, nmc, alpha, p, S, scale_matrix, G_mat_adj,
    matrix_accumulator_gibbs, start_point_first_gibbs
  )

  return(list(post_mean_omega=ans_hw[[1]], MC_average_Equation_9=ans_hw[[2]]))
}