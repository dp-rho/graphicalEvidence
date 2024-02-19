
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
  
  # Dimension of this iteration of Hao Wang
  p <- nrow(S)
  
  # Storage for sampled omega values
  omega_save <- array(0, dim = c(p, p, nmc))
  
  ind_noi_all = matrix(0, nrow=(p - 1), ncol=p)
  for (i in 1:p) {
    
    # Create ith col
    if (i == 1) {
      ind_noi <- 2:p
    }
    
    else if (i == p) {
      ind_noi <- 1:(p - 1)
    }
    
    else {
      ind_noi <- c(1:(i - 1), (i + 1):p)
    }
    
    ind_noi_all[, i] <- ind_noi
  }
  
  # These three initializations below are for computing the Normal density
  # in the evaluation of the term IV_{p-j+1}
  vec_log_normal_density <- rep(0, nmc)
  G_mat_last_col <- G_mat_adj[1:(p - 1), p]
  which_ones <- as.integer(G_mat_last_col == 1)
  logi_which_ones <- which(which_ones == 1)
  
  if (sum(which_ones) >= 1) {
    inv_C_required_store <- array(
      0, dim=c(sum(which_ones), sum(which_ones), nmc)
    )
    mean_vec_store <- matrix(0, nrow=sum(which_ones), ncol=nmc)
  }
  
  # Set initial value of Omega
  Omega <- start_point_first_gibbs
  
  
  # Start sampling
  for (iter in 1:(burnin + nmc)) {
    
    # Gibbs sampler with with Hao-Wang's decomposition
    for (i in 1:p) {
      V_mat_22 <- scale_matrix[i, i]
      ind_noi <- ind_noi_all[, i]
      V_mat_12 <- scale_matrix[ind_noi, i]
      
      S_21 <- S[ind_noi, i]
      S_22 <- S[i, i]
      vec_accumulator_21 <- -1 * matrix_accumulator_gibbs[ind_noi, i]
      
      # Sample gamma and beta
      gamma_sample <- rgamma(
        1, alpha + (n / 2) + 1, scale=(2 / (S_22 + V_mat_22))
      )
      inv_Omega_11 <- solve(Omega[ind_noi, ind_noi])
      inv_C <- (S_22 + V_mat_22) * inv_Omega_11

      # Identify which elements in beta to sample
      G_mat_current_col <- G_mat_adj[ind_noi, i]
      which_ones <- as.integer(G_mat_current_col == 1)
      logi_which_ones <- which(which_ones == 1)
      logi_which_zeros <- which(G_mat_current_col == 0)
      
      if (sum(which_ones) >= 1) {
        inv_C_required <- inv_C[logi_which_ones, logi_which_ones]
        inv_C_chol_required <- chol(inv_C_required)
        
        V_mat_12_required <- V_mat_12[logi_which_ones]
        S_21_required <- S_21[logi_which_ones]
        
        if (length(logi_which_zeros)) {
          
          vec_accumulator_21_required <- vec_accumulator_21[logi_which_zeros]
          inv_C_not_required <- inv_C[logi_which_zeros, logi_which_ones]
          
          vec_accumulator_21_required_mod <- t(
            t(vec_accumulator_21_required) %*% inv_C_not_required
          )

          mu_i_reduced <- -solve(
            inv_C_required,
            V_mat_12_required + S_21_required + vec_accumulator_21_required_mod
          )
          
          if ((iter > burnin) & (i == p)) {
            inv_C_required_store[, , iter - burnin] <- inv_C_required
            mean_vec_store[, iter - burnin] <- mu_i_reduced
          }
          
          beta_reduced <- mu_i_reduced + (
            solve(inv_C_chol_required, rnorm(sum(which_ones)))
          )
        }
        else {
          
          mu_i_reduced <- -solve(
            inv_C_required, V_mat_12_required + S_21_required
          )
          
          if ((iter > burnin) & (i == p)) {
            inv_C_required_store[, , iter - burnin] <- inv_C_required
            mean_vec_store[, iter - burnin] <- mu_i_reduced
          }
          
          beta_reduced <- mu_i_reduced + (
            solve(inv_C_chol_required, rnorm(sum(which_ones)))
          )
        }
        
        beta <- matrix(0, nrow=(p - 1), ncol=1)
        
        beta[logi_which_ones, 1] <- beta_reduced
        
        if (!(sum(which_ones) == (p - 1))) {
          beta[logi_which_zeros, 1] <- vec_accumulator_21_required
        }
      }
      else {
        beta <- vec_accumulator_21
      }

      omega_12 <- beta
      omega_22 <- gamma_sample + (t(beta) %*% inv_Omega_11 %*% beta)
      
      # Update omega
      Omega[i, ind_noi] <- omega_12
      Omega[ind_noi, i] <- omega_12
      Omega[i, i] <- omega_22
    }
    
    if (iter > burnin) {
      omega_save[, , iter - burnin] <- as.matrix(Omega)
    }
    
  }
  
  # Take accumulated mean of Omega from last nmc iterations
  post_mean_omega <- apply(omega_save, c(1, 2), mean)
  
  # Initialize time to calculate vec_log_normal_density
  cur_calc_time <- proc.time()
  
  ###################################
  ## RcppArmadillo implementation ###
  ###################################
  
  # RcppArmadillo implementation
  ans_ls <-calc_eq_9(
    which_ones, logi_which_ones, post_mean_omega, inv_C_required_store,
    mean_vec_store, p, nmc
  )

  ### R code time profiling ###
  calc_time <- proc.time() - cur_calc_time
  g_time_env$vec_log_normal_calc_time <- (
    g_time_env$vec_log_normal_calc_time + calc_time
  )
  ######################

  return(list(post_mean_omega=ans_ls[[1]], MC_average_Equation_9=ans_ls[[2]]))

  ###################################
  ## Native R implementation      ###
  ###################################
  
  # Calculate vec_log_normal_density by iterating through sampled values
  # if relevant, if no non zero values in G_mat_last_col, mean Eq9 is 0
  if (!sum(which_ones)){
    MC_average_Equation_9 <- 0
  }
  
  else {
    for (sample_index in 1:nmc) {
      inv_C_required <- inv_C_required_store[, , sample_index]
      mean_vec <- mean_vec_store[, sample_index]
      vec_log_normal_density[sample_index] <- log(
        mvtnorm::dmvnorm(
          post_mean_omega[p, logi_which_ones], mean_vec, solve(inv_C_required)
        )
      )
    }
    MC_average_Equation_9 <- log(mean(exp(vec_log_normal_density)))
  }
  
  ### R code time profiling ###
  calc_time <- proc.time() - cur_calc_time
  g_time_env$vec_log_normal_calc_time <- (
    g_time_env$vec_log_normal_calc_time + calc_time
  )
  ######################

  return(list(
    post_mean_omega=post_mean_omega,
    MC_average_Equation_9=MC_average_Equation_9
  ))
    
}