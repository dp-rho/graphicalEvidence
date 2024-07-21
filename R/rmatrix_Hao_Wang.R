
rmatrix_Hao_Wang <- function(
  S,
  n,
  burnin,
  nmc,
  prior,
  dof = 0,
  matrix_acc_gibbs = 0,
  lambda = 0
) {

  coded_prior <- switch(
    prior,
    'Wishart' = 0,
    'BGL' = 1,
    'GHS'= 2
  )

  p <- nrow(S)
  
  ##################################
  # RcppArmadillo implementation  ##
  ##################################

  hw_results <- mcmc_hw_rmatrix(
    n, burnin, nmc, p, coded_prior, dof, lambda, S, matrix_acc_gibbs
  )
  
  ##################################

  return(list(
    post_mean_omega=hw_results[[1]],
    post_mean_tau=hw_results[[2]],
    MC_avg_eq_9=hw_results[[3]]
  ))
}