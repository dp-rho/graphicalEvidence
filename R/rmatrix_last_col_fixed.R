
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
  
  # bind_random_samples_rmatrix(
  #   gamma_vec, rnorm_vec, rgig_vec, c(1)
  # )

  # RcppArmadillo implementation
  ans_last_col <- mcmc_last_col_rmatrix(
    n, burnin, nmc, p, dof, lambda, coded_prior, fixed_last_col, S,
    post_mean_tau, matrix_acc_gibbs, post_mean_omega
  )

  ##################################

  ### R code time profiling ###
  calc_time <- proc.time() - start_time_mcmc_last_col
  g_time_env$mcmc_last_col_calc_time <- (
    g_time_env$mcmc_last_col_calc_time + calc_time
  )
  ######################

  return(list(
    MC_avg_eq_11=ans_last_col[[1]],
    start_point_first_gibbs=ans_last_col[[2]],
    post_mean_omega_22=ans_last_col[[3]]
  ))

  ######################
  
  # Random generation debug
  # rnorm_index <- 1
  # rgig_index <- 1
  # rgamma_index <- 1
  # gamma_vec <- numeric((burnin + nmc) * p * (p - 1))
  # rnorm_vec <- numeric((burnin + nmc) * p * (p - 1))
  # rgig_vec <- numeric((burnin + nmc) * p * (p - 1))
  # runi_vec <- c()
  
  # omega_reduced_save <- array(0, dim=c(p, p, nmc))
  # gamma_subtractors <- numeric(nmc)
  # p_reduced <- p - 1
  # S_reduced <- as.matrix(S[1:p_reduced, 1:p_reduced])
  # matrix_acc_gibbs_reduced <- as.matrix(
  #   matrix_acc_gibbs[1:p_reduced, 1:p_reduced]
  # )
  # 
  # ind_noi_all = matrix(0, nrow=(p_reduced - 1), ncol=p_reduced)
  # if (p_reduced != 1) { 
  #   for (i in 1:p_reduced) {
  #     
  #     # Create ith col
  #     if (i == 1) {
  #       ind_noi <- 2:p_reduced
  #     }
  #     
  #     else if (i == p_reduced) {
  #       ind_noi <- 1:(p_reduced - 1)
  #     }
  #     
  #     else {
  #       ind_noi <- c(1:(i - 1), (i + 1):p_reduced)
  #     }
  #     ind_noi_all[,i] <- ind_noi
  #   }
  # }
  # 
  # omega_reduced <- as.matrix(post_mean_omega[1:p_reduced, 1:p_reduced])
  # tau_reduced <- as.matrix(post_mean_tau[1:p_reduced, 1:p_reduced])
  # if (prior == 'GHS') {
  #   nu_reduced <- matrix(1, nrow=p_reduced, ncol=p_reduced)
  # }
  # 
  # # Get prior specific parameters for gamma density
  # if (prior == 'Wishart') {
  #   shape_const <- (dof + n - p + 1) / 2
  #   scale_const <- 2 / (S[p, p] + 1)
  # }
  # else if (prior == 'BGL') {
  #   shape_const <- n / 2 + 1
  #   scale_const <- 2 / (S[p, p] + lambda)
  # }
  # else if (prior == 'GHS') {
  #   shape_const <- n / 2 + 1
  #   scale_const <- 2 / (S[p, p] + 1 / lambda)
  # }
  # 
  # # Begin main MCMC loop
  # for (iter in 1:(burnin + nmc)) {
  #   
  #   # First we update omega_pp which is nothing but sampling omega_22
  #   # with omega_12 held fixed
  #   inv_omega_11 <- solve(omega_reduced)
  # 
  #   # Sample omega_pp
  #   gamma_param <- rgamma(1, shape_const, 1 / scale_const)
  #   # gamma_vec[rgamma_index] <- gamma_param
  #   # rgamma_index <- rgamma_index + 1
  #   
  #   gamma_subtractor <- (
  #     t(fixed_last_col) %*% inv_omega_11 %*% fixed_last_col
  #   )[1]
  #   if (iter > burnin) {
  #     gamma_subtractors[iter - burnin] <- gamma_subtractor
  #   }
  #   omega_pp <- gamma_param + gamma_subtractor
  #   
  #   # cat("R CODE OMEGA PP:", omega_pp, "\n")
  #   
  #   # General case where we are not on the last iteration
  #   if (p_reduced != 1) {
  #     
  #     fixed_cross <- fixed_last_col %*% t(fixed_last_col)
  #     temp_matrix_acc <- (1 / omega_pp) * fixed_cross
  #     omega_reduced_tilde <- omega_reduced - temp_matrix_acc
  #     
  #     for (i in 1:p_reduced) {
  #       
  #       ind_noi <- ind_noi_all[, i]
  #       S_21_tilde <- S_reduced[ind_noi, i]
  #       
  #       if (prior == 'Wishart') {
  #         scale <- 2 / (S_reduced[i, i] + 1)
  #       }
  #       else if (prior == 'BGL') {
  #         vec_acc_21 <- matrix_acc_gibbs_reduced[ind_noi, i]
  #         temp_vec_acc_21 <- temp_matrix_acc[ind_noi, i]
  #         tau_12 <- tau_reduced[ind_noi, i]
  #         scale <- 2 / (S_reduced[i, i] + lambda)
  #       }
  #       else if (prior == 'GHS') {
  #         vec_acc_21 <- matrix_acc_gibbs_reduced[ind_noi, i]
  #         temp_vec_acc_21 <- temp_matrix_acc[ind_noi, i]
  #         tau_12 <- tau_reduced[ind_noi, i]
  #         nu_12 <- nu_reduced[ind_noi, i]
  #         scale <- 2 / (S_reduced[i, i] + 1 / lambda)
  #       }
  # 
  #       # Sampling from the gamma density of EQ 15
  #       gamma_param_tilde <- rgamma(1, shape_const, 1 / scale)
  #       # gamma_vec[rgamma_index] <- gamma_param_tilde
  #       # rgamma_index <- rgamma_index + 1
  #       # cat(i, "th gamma: ", gamma_param_tilde, '\n')
  #       
  #       tilde_w_11 <- omega_reduced_tilde[ind_noi, ind_noi]
  #       inv_omega_11 <- solve(tilde_w_11)
  #       # cat('inv omega 11 R: \n')
  #       # print(inv_omega_11)
  #       
  #       if (prior == 'Wishart') {
  #         inv_c <- (1 + S_reduced[i, i]) * inv_omega_11
  #         mu_i <- -solve(inv_c, S_21_tilde)
  #       }
  #       else if (prior == 'BGL') {
  #         inv_c <- (
  #           diag(1 / tau_12, nrow=(p_reduced - 1)) + 
  #           ((S_reduced[i, i] + lambda)) * inv_omega_11
  #         )
  #         mu_i <- -solve(
  #           inv_c, S_21_tilde + vec_acc_21 / tau_12 + temp_vec_acc_21 / tau_12
  #         )
  #       }
  #       else if (prior == 'GHS') {
  #         inv_c <- (
  #           diag(1 / (tau_12 * lambda^2), nrow=(p_reduced - 1)) + 
  #             (S_reduced[i, i] + 1 / lambda) * inv_omega_11
  #         )
  #         mu_i <- -solve(
  #           inv_c, (
  #             S_21_tilde + vec_acc_21 / (tau_12 * lambda^2) + 
  #             temp_vec_acc_21 / (tau_12 * lambda^2)
  #           )
  #         )
  #         # cat("R code solve for (?algebra): \n", S_21_tilde + vec_acc_21 / (tau_12 * lambda^2) + 
  #         #       temp_vec_acc_21 / (tau_12 * lambda^2), "\n")
  #       }
  #       
  #       cur_rnorm <- rnorm(p_reduced - 1)
  #       # rnorm_vec[rnorm_index:(rnorm_index + p_reduced - 2)] <- cur_rnorm
  #       # rnorm_index <- rnorm_index + p_reduced - 1
  #       
  #       beta <- mu_i + solve(chol(inv_c), cur_rnorm)
  #       # cat("R code beta: \n", beta, "\n")
  #       
  #       # Update omega tilde
  #       omega_reduced_tilde[i, ind_noi] <- beta
  #       omega_reduced_tilde[ind_noi, i] <- beta
  #       omega_reduced_tilde[i, i] <- gamma_param_tilde + (
  #         t(beta) %*% inv_omega_11 %*% beta
  #       )
  #       # cat("R code omega: \n")
  #       # print(omega_reduced_tilde)
  #       
  #       if (prior != 'Wishart') {
  #         mu_prime <- sqrt(lambda^2 / (beta + vec_acc_21 + temp_vec_acc_21)^2)
  #         
  #         # Update tau and nu if needed
  #         if (prior == 'BGL') {
  #         
  #           # Implementation of Generalized Inverse Gaussian sampler
  #           a_gig_tau <- lambda^2 / mu_prime^2
  #           u_12 <- numeric(p_reduced - 1)
  #           for (tau_id in 1:(p_reduced - 1)) {
  #             u_12[tau_id] <- gigrnd(-1 / 2, a_gig_tau[tau_id], lambda^2, 1)
  #           }
  #           
  #           tau_12 <- 1 / u_12
  #         }
  #         else if (prior == 'GHS') {
  #           rate <- beta^2 / (2 * lambda^2) + 1 / nu_12
  #           
  #           gammas1 <- rgamma(length(rate), 1, rate)
  #           # tau_12 <- 1 / rgamma_compiled(length(rate), 1, vec_rates=rate)
  #           tau_12 <- 1 / gammas1
  #           
  #           gammas2 <- rgamma(length(rate), 1, 1 + 1 / tau_12)
  #           # nu_12 <- 1 / rgamma_compiled(length(rate), 1, vec_rates=(1 + 1 / tau_12))
  #           nu_12 <- 1 / gammas2
  #           
  #           # last_index <- rgamma_index
  #           # for (j in 1:(p_reduced - 1)) {
  #           #   gamma_vec[rgamma_index + (2 * (j - 1))] <- gammas1[j]
  #           #   gamma_vec[rgamma_index + (2 * (j - 1)) + 1] <- gammas2[j]
  #           # }
  #           # rgamma_index <- rgamma_index + (2 * (p_reduced - 1))
  #           # cat(last_index, ":", rgamma_index)
  #           # print(gamma_vec[last_index:rgamma_index])
  #         }
  #         
  #         tau_reduced[i, ind_noi] <- tau_12
  #         tau_reduced[ind_noi, i] <- tau_12
  # 
  #         if (prior == 'GHS') {
  #           nu_reduced[i, ind_noi] <- nu_12
  #           nu_reduced[ind_noi, i] <- nu_12
  #         }
  #         # cat("R code tau: \n")
  #         # print(tau_reduced)
  #         # cat("R code nu: \n")
  #         # print(nu_reduced)
  #         # browser()
  #       }
  #     }
  #     
  #     omega_reduced <- omega_reduced_tilde + temp_matrix_acc
  #   }
  #   
  #   # Simple case for last iteration
  #   else {
  #     
  #     # Optimize by calculating once at the start
  #     if (prior == 'Wishart') {
  #       scale <- 2 / (S_reduced[1, 1] + 1)
  #     }
  #     else if (prior == 'BGL') {
  #       scale <- 2 / (S_reduced[1, 1] + lambda)
  #     }
  #     else if (prior == 'GHS') {
  #       scale <- 2 / (S_reduced[1, 1] + 1 / lambda)
  #     }
  #     gamma_param <- rgamma(1, shape_const, 1 / scale)
  #     omega_reduced[1, 1] <- gamma_param + (fixed_last_col[1]^2 / omega_pp)
  #   }
  #   
  #   if (iter > burnin) {
  #     omega_reduced_save[1:p_reduced, 1:p_reduced, iter - burnin] <- omega_reduced
  #     # OPTIMIZE CALCS: can assign fixed_last_col after full MCMC only once
  #     omega_reduced_save[p, 1:p_reduced, iter - burnin] <- fixed_last_col
  #     omega_reduced_save[p, 1:p_reduced, iter - burnin] <- fixed_last_col
  #     omega_reduced_save[p, p, iter - burnin] <- omega_pp
  #   }
  # }
  # 
  # post_mean_omega_22 <- mean(omega_reduced_save[p, p, ])
  # ind_noi <- 1:p_reduced
  # gamma_density_acc <- 0
  # 
  # # cat("R code post meant omega 22:", post_mean_omega_22, "\n")
  # for (sample in 1:nmc) {
  #   
  #   temp_gamma <- post_mean_omega_22 - gamma_subtractors[sample]
  #   
  #   # cat("R code gamma sub", sample ,": ", gamma_subtractors[sample], "\n")
  #   
  #   if (temp_gamma > 0) {
  #     gamma_density_acc <- gamma_density_acc + dgamma(
  #       temp_gamma, shape_const, scale=scale_const
  #     )
  #   }
  #   
  # }
  # 
  # # ### R code time profiling ###
  # calc_time <- proc.time() - start_time_mcmc_last_col
  # g_time_env$mcmc_last_col_calc_time <- (
  #   g_time_env$mcmc_last_col_calc_time + calc_time
  # )
  # cat("LAST COL CALC TIME: \n")
  # print(calc_time)
  # # ######################
  # 
  # return(list(
  #   post_mean_omega_22=post_mean_omega_22,
  #   MC_avg_eq_11=log(gamma_density_acc / nmc)
  # ))
  
}