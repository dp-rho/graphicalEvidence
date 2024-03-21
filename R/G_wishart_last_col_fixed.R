
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
  
  G_mat_adj_copy <- G_mat_adj
  matrix_accumulator_gibbs_copy <- matrix_accumulator_gibbs
  
  # gamma_subtractors <- numeric(nmc)
  # omega_reduced_acc <- matrix(0, nrow=(p - 1), ncol=(p - 1))
  # omega_22_acc <- 0
  # 
  # p_reduced <- p - 1
  # S_reduced <- as.matrix(S[1:p_reduced, 1:p_reduced])
  # 
  # matrix_accumulator_gibbs <- as.matrix(
  #   matrix_accumulator_gibbs[1:p_reduced, 1:p_reduced]
  # )
  # scale_matrix_required <- as.matrix(scale_matrix[1:p_reduced, 1:p_reduced])
  # G_mat_adj <- as.matrix(G_mat_adj[1:p_reduced, 1:p_reduced])
  # 
  # if (p_reduced != 1) {
  #   ind_noi_all <- matrix(0, nrow=(p_reduced - 1), ncol=p_reduced)
  #   for (i in 1:p_reduced) {
  #     if (i == 1) {
  #       ind_noi <- 2:p_reduced
  #     }
  #     else if (i == p_reduced) {
  #       ind_noi <- 1:(p_reduced - 1)
  #     }
  #     else {
  #       ind_noi <- c(1:(i - 1), (i + 1):p_reduced)
  #     }
  #     ind_noi_all[, i] <- ind_noi
  #   }
  # }
  # 
  # omega_reduced <- post_mean_omega[1:p_reduced, 1:p_reduced]
  # 
  # # Save random sampled values to export to C++ for validation
  # # gamma_vals <- numeric((burnin + nmc) * p)
  # # gamma_index <- 1
  # # norm_vals <- c((burnin + nmc) * p^2)
  # # norm_index <- 1
  # 
  # # Start MCMC sampling
  # for (iter in 1:(burnin + nmc)) {
  # 
  #   # First update omega_pp
  #   # inv_omega_11 <- solve(omega_reduced)
  #   s_22 <- S[p, p]
  #   V_mat_22 <- scale_matrix[p, p]
  # 
  #   # Sample omega_22
  #   gamma_sample <- rgamma(1, alpha + (n / 2) + 1, scale=(2 / (s_22 + V_mat_22)))
  #   omega_pp <- gamma_sample + (
  #     t(fixed_last_col) %*% solve(omega_reduced, fixed_last_col)
  #   )[1]
  # 
  #   # Saving random samples for replication in compiled code
  #   # gamma_vals[gamma_index] <- gamma_sample
  #   # gamma_index <- gamma_index + 1
  #   #####
  # 
  #   if (p_reduced != 1) {
  #     omega_reduced_tilde <- omega_reduced + (
  #       (1 / omega_pp) * (fixed_last_col %*% t(fixed_last_col))
  #     )
  # 
  #     temp_matrix_accumulator <- (1 / omega_pp) * (
  #       (fixed_last_col %*% t(fixed_last_col))
  #     )
  # 
  #     # sample omega_11_tilde
  #     for (i in 1:p_reduced) {
  #       V_mat_22 <- scale_matrix_required[i, i]
  #       ind_noi <- ind_noi_all[, i]
  #       V_mat_12 <- scale_matrix_required[ind_noi, i]
  # 
  #       s_21_tilde <- S_reduced[ind_noi, i]
  #       s_22_tilde <- S_reduced[i, i]
  # 
  #       # -1 is critical for crucial indicator function
  #       vec_accumulator_21 <- -1 * matrix_accumulator_gibbs[ind_noi, i]
  #       temp_vec_accumulator_21 <- -1 * temp_matrix_accumulator[ind_noi, i]
  # 
  #       gamma_sample <- rgamma(
  #         1, alpha + (n / 2) + 1, scale=(2 / (s_22_tilde + V_mat_22))
  #       )
  #       # Saving random samples for replication in compiled code
  #       # gamma_vals[gamma_index] <- gamma_sample
  #       # gamma_index <- gamma_index + 1
  #       #####
  # 
  #       tilde_w_11 <- omega_reduced_tilde[ind_noi, ind_noi]
  #       inv_omega_11 <- solve(tilde_w_11)
  # 
  #       inv_C <-(s_22_tilde * V_mat_22) * inv_omega_11
  # 
  #       G_mat_cur_col <- G_mat_adj[ind_noi, i]
  #       which_ones <- as.integer(G_mat_cur_col == 1)
  #       logi_which_ones <- which(which_ones == 1)
  #       which_zeros <- as.integer(G_mat_cur_col == 0)
  #       logi_which_zeros <- which(which_zeros == 1)
  # 
  #       if (sum(which_ones) >= 1) {
  #         inv_C_required <- inv_C[logi_which_ones, logi_which_ones]
  #         inv_C_chol_required <- chol(inv_C_required)
  # 
  #         V_mat_12_required <- V_mat_12[logi_which_ones]
  #         s_21_tilde_required <- s_21_tilde[logi_which_ones]
  # 
  #         if (sum(which_zeros) >= 1) {
  #           vec_accumulator_21_required <- vec_accumulator_21[logi_which_zeros]
  #           inv_C_not_required <- inv_C[logi_which_zeros, logi_which_ones]
  # 
  #           vec_accumulator_21_required_mod <- as.vector(
  #             t(vec_accumulator_21_required) %*% inv_C_not_required
  #           )
  # 
  #           temp_vec_accumulator_21_required <- (
  #             temp_vec_accumulator_21[logi_which_zeros]
  #           )
  # 
  #           temp_vec_accumulator_21_required_mod <- as.vector(
  #             temp_vec_accumulator_21_required %*% inv_C_not_required
  #           )
  # 
  #           # cat(paste0("i: ", i))
  #           # print("inv c required: ")
  #           # print(inv_C_required)
  #           # print("solve for: ")
  #           # print(V_mat_12_required + s_21_tilde_required +
  #           #         vec_accumulator_21_required_mod +
  #           #         temp_vec_accumulator_21_required_mod)
  # 
  #           mu_i_reduced <- -solve(
  #             inv_C_required,
  #             V_mat_12_required + s_21_tilde_required +
  #             vec_accumulator_21_required_mod +
  #             temp_vec_accumulator_21_required_mod
  #           )
  #         }
  #         else {
  #           rnorm_sample <- rnorm(sum(which_ones))
  #           # Saving random samples for replication in compiled code
  #           # norm_vals[norm_index:(norm_index + length(rnorm_sample) - 1)] <- rnorm_sample
  #           # norm_index <- norm_index + length(rnorm_sample)
  #           ####
  #           mu_i_reduced <- -solve(inv_C_required, rnorm_sample)
  #         }
  # 
  #         rnorm_sample <- rnorm(sum(which_ones))
  #         # Saving random samples for replication in compiled code
  #         # norm_vals[norm_index:(norm_index + length(rnorm_sample) - 1)] <- rnorm_sample
  #         # norm_index <- norm_index + length(rnorm_sample)
  #         ####
  # 
  #         beta_reduced <- mu_i_reduced + (
  #           solve(inv_C_chol_required, rnorm_sample)
  #         )
  # 
  #         beta <- rep(0, p_reduced - 1)
  #         beta[logi_which_ones] <- beta_reduced
  # 
  #         if (sum(which_ones) != (p_reduced - 1)) {
  #           beta[logi_which_zeros] <- (
  #             vec_accumulator_21_required + temp_vec_accumulator_21_required
  #           )
  #         }
  #       }
  #       else {
  #         beta <- vec_accumulator_21 + temp_vec_accumulator_21
  #       }
  # 
  #       omega_12 <- beta
  #       omega_22 <- gamma_sample + (
  #         t(beta) %*% inv_omega_11 %*% beta
  #       )
  # 
  #       omega_reduced_tilde[i, ind_noi] <- omega_12
  #       omega_reduced_tilde[ind_noi, i] <- omega_12
  #       omega_reduced_tilde[i, i] <- omega_22
  #     }
  # 
  #     omega_reduced <- omega_reduced_tilde + (
  #       (1 / omega_pp) * (fixed_last_col %*% t(fixed_last_col))
  #     )
  #   }
  #   else {
  #     s_22 <- S_reduced[1, 1]
  #     V_mat_22 <- scale_matrix_required[1, 1]
  #     gamma_sample <- rgamma(1, alpha + (n / 2) + 1, scale=(2 / (s_22 + V_mat_22)))
  #     # Saving random samples for replication in compiled code
  #     # gamma_vals[gamma_index] <- gamma_sample
  #     # gamma_index <- gamma_index + 1
  #     #####
  # 
  #     # print("ELSE")
  #     # print(fixed_last_col)
  #     # print(omega_reduced)
  # 
  #     omega_reduced <- gamma_sample + matrix(
  #       fixed_last_col[1]^2 / omega_pp, nrow=1, ncol=1
  #     )
  #     # omega_reduced <- gamma_sample + (
  #     #   t(fixed_last_col) %*% solve(omega_pp, fixed_last_col)
  #     # )
  #     # print("GAMMA SAMPLE")
  #     # print(gamma_sample)
  #     # print("OMEGA REDUCED POST")
  #     # print(omega_reduced)
  # 
  #   }
  # 
  #   if (iter > burnin) {
  #     gamma_subtractors[iter - burnin] <- (
  #       # t(fixed_last_col) %*% solve(omega_reduced) %*% fixed_last_col
  #       t(fixed_last_col) %*% solve(omega_reduced, fixed_last_col)
  #     )
  #     omega_22_acc <- omega_22_acc + omega_pp
  #     omega_reduced_acc <- omega_reduced_acc + as.matrix(omega_reduced)
  #   }
  # }
  # 
  # omega_22_acc <- omega_22_acc / nmc
  # start_point_first_gibbs <- omega_reduced_acc / nmc
  
  # print("R ANS")
  # print(gamma_subtractors)
  # print(start_point_first_gibbs)
  # print(omega_22_acc)
  
  # vec_gamma_density <- rep(-Inf, nmc)
  
  # bind_random_samples(gamma_vals, norm_vals)
  ##################################
  # RcppArmadillo implementation ###
  ##################################
  
  # RcppArmadillo implementation
  ans_last_col <- mcmc_last_col(
    n, burnin, nmc, alpha, p, fixed_last_col, S, scale_matrix, G_mat_adj_copy,
    matrix_accumulator_gibbs_copy, post_mean_omega
  )
  
  # print(ans_last_col)
  
  compiled_gamma_subtractors <- ans_last_col[[1]]
  compiled_omega_reduced_acc <- ans_last_col[[2]]
  compiled_omega_22_acc <- ans_last_col[[3]]
  
  ##################################
  

  # Initialize time to calculate vec_log_normal_density
  start_time_vec <- proc.time()
  
  ###################################
  ## RcppArmadillo implementation ###
  ###################################
  
  # RcppArmadillo implementation
  MC_average_Equation_11 <- calc_eq_11(
    compiled_omega_22_acc, S[p, p], scale_matrix[p, p], alpha, # compiled_omega_22_acc
    n, nmc, compiled_gamma_subtractors # compiled_gamma_subtractors
  )

  
  ###################################
  ## Native R implementation      ###
  ###################################

  # for (sample_index in 1:nmc) {
  # 
  #   temp_gamma <- omega_22_acc - gamma_subtractors[sample_index]
  # 
  #   if (temp_gamma > 0) {
  #     vec_gamma_density[sample_index] <- dgamma(
  #       temp_gamma, alpha + (n / 2) + 1, 
  #       scale=(2 / (S[p, p] + scale_matrix[p, p]))
  #     )
  #   }
  # }
  # 
  # MC_average_Equation_11 <- log(mean(vec_gamma_density))
  
  ###################################
  
  ### R code time profiling ###
  calc_time <- proc.time() - start_time_vec
  g_time_env$vec_gamma_calc_time <- (
    g_time_env$vec_gamma_calc_time + calc_time
  )
  calc_time <- proc.time() - start_time_mcmc_last_col
  g_time_env$mcmc_last_col_calc_time <- (
    g_time_env$mcmc_last_col_calc_time + calc_time
  )
  ######################
  
  return(list(
    MC_average_Equation_11=MC_average_Equation_11,
    # start_point_first_gibbs=start_point_first_gibbs,
    start_point_first_gibbs=compiled_omega_reduced_acc,
    # post_mean_omega_22=omega_22_acc
    post_mean_omega_22=compiled_omega_22_acc
  ))
}