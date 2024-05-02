
BGL_Hao_Wang <- function(
  S,
  n,
  burnin,
  nmc,
  lambda,
  matrix_acc_gibbs
) {
  
  # Initialize storage
  p <- nrow(S)
  tau_acc <- 0
  omega_acc <- matrix(0, nrow=p, ncol=p)
  normal_density_acc <- 0
  inv_c_store <- array(0, dim=c(p - 1, p - 1, nmc))
  mean_vec_store <- matrix(0, nrow=(p - 1), ncol=nmc)
  
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
  
  omega <- diag(p)
  cur_sigma <- solve(omega)
  tau <- matrix(1, nrow=p, ncol=p) 
  
  # Begin main MCMC loop
  for (iter in 1:(burnin + nmc)) {
    
    # Gibb's sampler for Omega with Hao-Wang's decomposition
    for (i in 1:p) {
      ind_noi <- ind_noi_all[,i]
      sigma_11 <- cur_sigma[ind_noi, ind_noi]
      sigma_12 <- cur_sigma[ind_noi, i]
      sigma_22 <- cur_sigma[i, i]
      
      S_21 <- S[ind_noi, i]
      S_22 <- S[i, i]
      vec_acc_21 <- matrix_acc_gibbs[ind_noi, i]
      tau_12 <- tau[ind_noi, i]
      
      # Sampling from the gamma density of EQ 15
      gamma_param <- rgamma(1, n / 2 + 1, scale=2 / (S_22 + lambda))
      
      inv_omega_11 <- sigma_11 - sigma_12 %*% t(sigma_12) / sigma_22
      
      inv_c <- diag(1 / tau_12, nrow=(p - 1)) + ((S_22 + lambda) * inv_omega_11)
      mu_i <- -solve(
        inv_c, (S_21 + (vec_acc_21 / tau_12))
      )
      
      # If burnin is complete, save 
      if ((iter > burnin) & (i == p)) {
        inv_c_store[, , (iter - burnin)] <- inv_c
        mean_vec_store[, (iter - burnin)] <- mu_i
      }
      
      # Sampling from the Normal density of Equation (15) in the paper
      beta <- mu_i + solve(chol(inv_c), rnorm(p - 1))
      
      omega_22 <- gamma_param + t(beta) %*% inv_omega_11 %*% beta
      
      # Update omega
      omega[i, ind_noi] <- beta
      omega[ind_noi, i] <- beta
      omega[i, i] <- omega_22
      
      # TODO
      omega_beta <- inv_omega_11 %*% beta;
      sigma_11 <- inv_omega_11 + (omega_beta %*% t(omega_beta) / gamma_param)
      sigma_12 <- -omega_beta / gamma_param
      sigma_22 <- 1 / gamma_param
      cur_sigma[ind_noi, ind_noi] <- sigma_11
      cur_sigma[ind_noi, i] <- sigma_12
      cur_sigma[i, ind_noi] <- sigma_12
      cur_sigma[i, i] <- sigma_22
      
      # Update tau
      mu_prime_sq <- lambda^2 / (beta + vec_acc_21)^2
      a_gig_tau <- lambda^2 / mu_prime_sq
      u_12 <- numeric(p - 1)
      for (tau_id in 1:(p - 1)){
        u_12[tau_id] <- gigrnd(-1 / 2, a_gig_tau[tau_id], lambda^2, 1)
      }
      tau_12 <- 1 / u_12
      
      tau[i, ind_noi] <- tau_12
      tau[ind_noi, i] <- tau_12
    }
    
    if (iter > burnin) {
      omega_acc <- omega_acc + omega
      tau_acc <- tau_acc + tau
    }
  }
  
  ind_noi <- ind_noi_all[, p]
  omega_acc <- omega_acc / nmc
  tau_acc <- tau_acc / nmc

  for (sample_index in 1:nmc) {
    normal_density_acc <- normal_density_acc + mvtnorm::dmvnorm(
      omega_acc[p, ind_noi], mean_vec_store[, sample_index], 
      solve(inv_c_store[, , sample_index])
    )
  }
  MC_avg_eq_9 <- log(normal_density_acc / nmc)
  return(list(
    MC_avg_eq_9=log(normal_density_acc / nmc),
    post_mean_omega=omega_acc,
    post_mean_tau=tau_acc
  ))
}