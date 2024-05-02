
# Last col restricted sampler for BGL
BGL_last_col_fixed <- function(
  S,
  n,
  burnin,
  nmc,
  lambda,
  fixed_last_col,
  matrix_acc_gibbs,
  post_mean_omega,
  post_mean_tau
) {
  
  # Initialize storage
  p <- nrow(S)
  omega_reduced_save <- array(0, dim=c(p, p, nmc))
  gamm_val_store <- numeric(nmc)
  p_reduced <- p - 1
  S_reduced <- as.matrix(S[1:p_reduced, 1:p_reduced])
  matrix_acc_gibbs_reduced <- as.matrix(
    matrix_acc_gibbs[1:p_reduced, 1:p_reduced]
  )
  
  ind_noi_all = matrix(0, nrow=(p_reduced - 1), ncol=p_reduced)
  if (p_reduced != 1) { 
    for (i in 1:p_reduced) {
      
      # Create ith col
      if (i == 1) {
        ind_noi <- 2:p_reduced
      }
      
      else if (i == p_reduced) {
        ind_noi <- 1:(p_reduced - 1)
      }
      
      else {
        ind_noi <- c(1:(i - 1), (i + 1):p_reduced)
      }
      ind_noi_all[,i] <- ind_noi
    }
  }
  
  omega_reduced <- as.matrix(post_mean_omega[1:p_reduced, 1:p_reduced])
  tau_reduced <- as.matrix(post_mean_tau[1:p_reduced, 1:p_reduced])
  
  # Begin main MCMC loop
  for (iter in 1:(burnin + nmc)) {
    
    # First we update omega_pp which is nothing but sampling omega_22
    # with omega_12 held fixed
    inv_omega_11 <- solve(omega_reduced)
    S_22 <- S[p, p]
    
    # Sample omega_pp
    gamma_param <- rgamma(1, n / 2 + 1, scale=2 / (S_22 + lambda))
    omega_pp <- gamma_param + (
      t(fixed_last_col) %*% inv_omega_11 %*% fixed_last_col
    )[1]
    
    # General case where we are not on the last iteration
    if (p_reduced != 1) {
      
      fixed_cross <- fixed_last_col %*% t(fixed_last_col)
      temp_matrix_acc <- (1 / omega_pp) * fixed_cross
      omega_reduced_tilde <- omega_reduced - temp_matrix_acc
      
      for (i in 1:p_reduced) {
        
        ind_noi <- ind_noi_all[, i]
        S_21_tilde <- S_reduced[ind_noi, i]
        S_22_tilde <- S_reduced[i, i]
        vec_acc_21 <- matrix_acc_gibbs_reduced[ind_noi, i]
        temp_vec_acc_21 <- temp_matrix_acc[ind_noi, i]
        tau_12 <- tau_reduced[ind_noi, i]
        # CONDENSE WITH ELSE
        gamma_param_tilde <- rgamma(1, n / 2 + 1, scale=2 / (S_22_tilde + lambda))
        
        tilde_w_11 <- omega_reduced_tilde[ind_noi, ind_noi]
        inv_omega_11 <- solve(tilde_w_11)
        inv_c <- (diag(1 / tau_12, nrow=(p_reduced - 1)) + ((S_22_tilde + lambda)) * inv_omega_11)
        
        mu_i <- -solve(
          inv_c, S_21_tilde + vec_acc_21 / tau_12 + temp_vec_acc_21 / tau_12
        )
        inv_c_chol <- chol(inv_c)
        beta <- mu_i + solve(inv_c_chol, rnorm(p_reduced - 1))
        
        # Update omega tilde
        omega_reduced_tilde[i, ind_noi] <- beta
        omega_reduced_tilde[ind_noi, i] <- beta
        omega_reduced_tilde[i, i] <- gamma_param_tilde + (
          t(beta) %*% inv_omega_11 %*% beta
        )
        
        mu_prime <- sqrt(lambda^2 / (beta + vec_acc_21 + temp_vec_acc_21)^2)
        
        # Implementation of Generalized Inverse Gaussian sampler
        a_gig_tau <- lambda^2 / mu_prime^2
        u_12 <- numeric(p_reduced - 1)
        for (tau_id in 1:(p_reduced - 1)) {
          u_12[tau_id] <- gigrnd(-1 / 2, a_gig_tau[tau_id], lambda^2, 1)
        }
        
        tau_12 <- 1 / u_12
        tau_reduced[i, ind_noi] <- tau_12
        tau_reduced[ind_noi, i] <- tau_12
      }
      
      omega_reduced <- omega_reduced_tilde + temp_matrix_acc
    }
    
    # Simple case for last iteration
    else {
      gamma_param <- rgamma(1, n / 2 + 1, scale=2 / (S_reduced[1, 1] + lambda))
      omega_reduced[1, 1] <- gamma_param + (fixed_last_col[1] ^ 2 / omega_pp)
    }
    
    if (iter > burnin) {
      omega_reduced_save[1:p_reduced, 1:p_reduced, iter - burnin] <- omega_reduced
      # OPTIMIZE CALCS: can assign fixed_last_col after full MCMC only once
      omega_reduced_save[p, 1:p_reduced, iter - burnin] <- fixed_last_col
      omega_reduced_save[p, 1:p_reduced, iter - burnin] <- fixed_last_col
      omega_reduced_save[p, p, iter - burnin] <- omega_pp
    }
  }
  
  post_mean_omega_22 <- mean(omega_reduced_save[p, p, ])
  ind_noi <- 1:p_reduced
  gamma_density_acc <- 0
  
  for (sample in 1:nmc) {
    
    temp_gamma <- post_mean_omega_22 - (
      t(fixed_last_col) %*% 
      solve(omega_reduced_save[ind_noi, ind_noi, sample]) %*%
      fixed_last_col
    )
    
    if (temp_gamma > 0) {
      gamma_density_acc <- gamma_density_acc + dgamma(
        temp_gamma, (n / 2) + 1, scale=(2 / (S[p, p] + lambda))
      )
    }
    else {
      browser()
      gamma_density_acc <- -Inf
      break
    }
  }
  
  return(list(
    post_mean_omega_22=post_mean_omega_22,
    MC_avg_eq_11=log(gamma_density_acc / nmc)
  ))
  
}