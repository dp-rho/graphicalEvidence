
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
  
  p <- nrow(S)
  omega_save <- array(0, dim=c(p, p, nmc))
  
  p_reduced <- p - 1 
  S_reduced <- as.matrix(S[1:p_reduced, 1:p_reduced])
  
  matrix_accumulator_gibbs <- as.matrix(
    matrix_accumulator_gibbs[1:p_reduced, 1:p_reduced]
  )
  scale_matrix_required <- as.matrix(scale_matrix[1:p_reduced, 1:p_reduced])
  G_mat_adj <- as.matrix(G_mat_adj[1:p_reduced, 1:p_reduced])
  
  if (p_reduced != 1) {
    ind_noi_all <- matrix(0, nrow=(p_reduced - 1), ncol=p_reduced)
    for (i in 1:p_reduced) {
      if (i == 1) {
        ind_noi <- 2:p_reduced
      }
      else if (i == p_reduced) {
        ind_noi <- 1:(p_reduced - 1)
      }
      else {
        ind_noi <- c(1:(i - 1), (i + 1):p_reduced)
      }
      ind_noi_all[, i] <- ind_noi
    }
  }
  
  omega_reduced <- post_mean_omega[1:p_reduced, 1:p_reduced]
  
  # Start MCMC sampling
  for (iter in 1:(burnin + nmc)) {
    
    # First update omega_pp
    inv_omega_11 <- solve(omega_reduced)
    s_22 <- S[p, p]
    V_mat_22 <- scale_matrix[p, p]
    
    # Sample omega_22
    gamma_sample <- rgamma(1, alpha + (n / 2) + 1, scale=(2 / (s_22 + V_mat_22)))
    omega_pp <- gamma_sample + (
      t(fixed_last_col) %*% inv_omega_11 %*% fixed_last_col
    )[1]

    if (p_reduced != 1) {
      omega_reduced_tilde <- omega_reduced + (
        (1 / omega_pp) * (fixed_last_col %*% t(fixed_last_col))
      )
      
      temp_matrix_accumulator <- (1 / omega_pp) * (
        (fixed_last_col %*% t(fixed_last_col))
      )
      
      # sample omega_11_tilde
      for (i in 1:p_reduced) {
        V_mat_22 <- scale_matrix_required[i, i]
        ind_noi <- ind_noi_all[, i]
        V_mat_12 <- scale_matrix_required[ind_noi, i]
        
        s_21_tilde <- S_reduced[ind_noi, i]
        s_22_tilde <- S_reduced[i, i]
        
        # -1 is critical for crucial indicator function
        vec_accumulator_21 <- -1 * matrix_accumulator_gibbs[ind_noi, i]
        temp_vec_accumulator_21 <- -1 * temp_matrix_accumulator[ind_noi, i]
        
        gamma_sample <- rgamma(
          1, alpha + (n / 2) + 1, scale=(2 / (s_22_tilde + V_mat_22))
        )

        tilde_w_11 <- omega_reduced_tilde[ind_noi, ind_noi]
        inv_omega_11 <- solve(tilde_w_11)
        
        inv_C <-(s_22_tilde * V_mat_22) * inv_omega_11
        
        G_mat_cur_col <- G_mat_adj[ind_noi, i]
        which_ones <- as.integer(G_mat_cur_col == 1)
        logi_which_ones <- which(which_ones == 1)
        which_zeros <- as.integer(G_mat_cur_col == 0)
        logi_which_zeros <- which(which_zeros == 1)
        
        if (sum(which_ones) >= 1) {
          inv_C_required <- inv_C[logi_which_ones, logi_which_ones]
          inv_C_chol_required <- chol(inv_C_required)
          
          V_mat_12_required <- V_mat_12[logi_which_ones]
          s_21_tilde_required <- s_21_tilde[logi_which_ones]
          
          if (sum(which_zeros) >= 1) {
            vec_accumulator_21_required <- vec_accumulator_21[logi_which_zeros]
            inv_C_not_required <- inv_C[logi_which_zeros, logi_which_ones]
            
            vec_accumulator_21_required_mod <- as.vector(
              t(vec_accumulator_21_required) %*% inv_C_not_required
            )
            
            temp_vec_accumulator_21_required <- (
              temp_vec_accumulator_21[logi_which_zeros]
            )
            
            temp_vec_accumulator_21_required_mod <- as.vector(
              temp_vec_accumulator_21_required %*% inv_C_not_required
            )

            mu_i_reduced <- -solve(
              inv_C_required,
              V_mat_12_required + s_21_tilde_required + 
              vec_accumulator_21_required_mod + 
              temp_vec_accumulator_21_required_mod
            )
            
            beta_reduced <- mu_i_reduced + (
              solve(inv_C_chol_required, rnorm(sum(which_ones)))
            )
          }
          else {
            mu_i_reduced <- -solve(inv_C_required, rnorm(sum(which_ones)))
            beta_reduced <- mu_i_reduced + (
              solve(inv_C_chol_required, rnorm(sum(which_ones)))
            )
          }
          
          beta <- rep(0, p_reduced - 1)
          beta[logi_which_ones] <- beta_reduced
          
          if (sum(which_ones) != (p_reduced - 1)) {
            beta[logi_which_zeros] <- (
              vec_accumulator_21_required + temp_vec_accumulator_21_required
            )
          }
        }
        else {
          beta <- vec_accumulator_21 + temp_vec_accumulator_21
        }
        
        omega_12 <- beta
        omega_22 <- gamma_sample + (
          t(beta) %*% inv_omega_11 %*% beta
        )
        
        omega_reduced_tilde[i, ind_noi] <- omega_12
        omega_reduced_tilde[ind_noi, i] <- omega_12
        omega_reduced_tilde[i, i] <- omega_22
      }
      
      omega_reduced <- omega_reduced_tilde + (
        (1 / omega_pp) * (fixed_last_col %*% t(fixed_last_col))
      )
    }
    else {
      s_22 <- S_reduced[1, 1]
      V_mat_22 <- scale_matrix_required[1, 1]
      gamma_sample <- rgamma(1, alpha + (n / 2) + 1, scale=(2 / (s_22 + V_mat_22)))
      omega_reduced <- gamma_sample + (
        t(fixed_last_col) %*% solve(omega_pp) %*% fixed_last_col
      )
      
    }
    
    if (iter < burnin) {
      omega_save[1:p_reduced, 1:p_reduced, iter - burnin] <- (
        as.matrix(omega_reduced)
      )
      omega_save[p, 1:p_reduced, iter - burnin] <- fixed_last_col
      omega_save[1:p_reduced, p , iter - burnin] <- fixed_last_col
      omega_save[p, p, iter - burnin] <- omega_pp
    }
  }
  
  post_mean_omega_22 <- mean(omega_save[p, p, ])
  vec_log_gamma_density <- rep(-Inf, nmc)
  ind_noi <- 1:p_reduced

  for (sample_index in 1:nmc) {
    inv_omega_11 <- solve(omega_save[ind_noi, ind_noi, sample_index])
    temp_gamma <- post_mean_omega_22 - (
      t(fixed_last_col) %*% inv_omega_11 %*% fixed_last_col
    )

    if (temp_gamma > 0) {
      vec_log_gamma_density[sample_index] <- log(dgamma(
        temp_gamma, alpha + (n / 2) + 1, 
        scale=(2 / (S[p, p] + scale_matrix[p, p]))
      ))
    }
  }
  
  MC_average_Equation_11 <- log(mean(exp(vec_log_gamma_density)))
  start_point_first_gibbs <- apply(omega_save, c(1, 2), mean)
  start_point_first_gibbs <- start_point_first_gibbs[1:p_reduced, 1:p_reduced]
  
  return(list(
    MC_average_Equation_11=MC_average_Equation_11,
    start_point_first_gibbs=start_point_first_gibbs,
    post_mean_omega_22=post_mean_omega_22
  ))
}