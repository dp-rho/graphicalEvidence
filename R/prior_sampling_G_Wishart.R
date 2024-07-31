
# Gets the initial starting point for gibbs sampler
prior_sampling_G_Wishart <- function(
    p, 
    burnin_prior, 
    nmc_prior, 
    G_mat_adj, 
    scale_matrix, 
    alpha
) {
  
  omega_save <- array(0, dim=c(p, p, nmc_prior))
  
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
    ind_noi_all[,i] <- ind_noi
  }
  
  # DEBUG
  vec_gammas <- numeric(p * (burnin_prior + nmc_prior))
  vec_gamma_index <- 1
  vec_norms <- numeric(p * p * (burnin_prior + nmc_prior))
  vec_norm_index <- 1
  
  # Set intial values
  Omega <- G_mat_adj
  Omega[Omega == 1] <- 0.01
  diag(Omega) <- 1
  for (iter in 1:(burnin_prior + nmc_prior)) {
    
    # Gibb's sampler for Omega with Hao-Wang's decomposition
    for (i in 1:p) {
      V_mat_22 <- scale_matrix[i, i]
      ind_noi <- ind_noi_all[, i]
      V_mat_12 <- scale_matrix[ind_noi, i]
      
      # Sample gamma and beta
      gamma_param <- rgamma(1, alpha + 1, scale=(2 / V_mat_22))
      
      # vec_gammas[vec_gamma_index] <- gamma_param
      # vec_gamma_index <- vec_gamma_index + 1
      
      inv_Omega_11 <- solve(Omega[ind_noi, ind_noi])
      inv_C <- V_mat_22 * inv_Omega_11
      
      # Identify which elements in beta to sample
      G_mat_current_col <- G_mat_adj[ind_noi, i]
      which_ones <- as.integer(G_mat_current_col == 1)
      logi_which_ones <- which(which_ones == 1)
      
      beta <- matrix(0, nrow=(p - 1), ncol=1)
      if (sum(which_ones) >= 1) {
        
        inv_C_required <- inv_C[logi_which_ones, logi_which_ones]
        inv_C_chol_required <- chol(inv_C_required)
        V_mat_12_required = V_mat_12[logi_which_ones]
        mu_i_reduced <- -solve(inv_C_required, V_mat_12_required)
        rnorm_sample <- rnorm(sum(which_ones))
        
        # vec_norms[vec_norm_index:(vec_norm_index + length(rnorm_sample) - 1)] <- rnorm_sample
        # vec_norm_index <- vec_norm_index + length(rnorm_sample)
        
        # print("mu i reduced target")
        # print(V_mat_12_required)
        # print("mu i reduced")
        # print(mu_i_reduced)
        # print("solve chol system")
        # print(solve(inv_C_chol_required, rnorm_sample))
        
        beta_reduced <- (
          mu_i_reduced + (
            solve(inv_C_chol_required, rnorm_sample)
          )
        )

        beta[logi_which_ones, 1] <- beta_reduced;
      }
      
      # cat(paste0(i, "th beta: \n"))
      # print(beta)
      
      omega_12 = beta; 
      omega_22 = gamma_param + t(beta) %*% inv_Omega_11 %*% beta;
      
      # cat("omega beta:\n")
      # print(inv_Omega_11 %*% beta)
      # cat("gamma param:", gamma_param, "\n")
      # print(inv_Omega_11)
      # cat("omega_22:", omega_22, "\n")
        
      # update Omega and Sigma
      Omega[i, ind_noi] <- omega_12
      Omega[ind_noi, i] <- omega_12
      Omega[i, i] = omega_22
    }
    
    # Save results past burnin
    if (iter > burnin_prior) {
      omega_save[, , (iter - burnin_prior)] <- as.matrix(Omega)
    }
  }
  # bind_random_samples(vec_gammas, vec_norms)
  
  return(omega_save)
}
