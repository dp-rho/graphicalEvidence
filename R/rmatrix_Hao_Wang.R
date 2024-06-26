
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
  
  # Initialize time to calculate Hao Wang sampler
  start_time_mcmc_hw <- proc.time()

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
  
  # bind_random_samples_rmatrix(
  #   gamma_vec, rnorm_vec, rgig_vec, c(1)
  # )

  hw_results <- mcmc_hw_rmatrix(
    n, burnin, nmc, p, coded_prior, dof, lambda, S, matrix_acc_gibbs
  )

  ### R code time profiling ###
  calc_time <- proc.time() - start_time_mcmc_hw
  g_time_env$mcmc_hw_calc_time <- (
    g_time_env$mcmc_hw_calc_time + calc_time
  )
  ######################

  return(list(
    post_mean_omega=hw_results[[1]],
    post_mean_tau=hw_results[[2]],
    MC_avg_eq_9=hw_results[[3]]
  ))
  
  ##################################
  
  # Initialize storage
  tau_acc <- matrix(0, nrow=p, ncol=p)
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
  if (prior == 'GHS'){
    nu <- matrix(1, nrow=p, ncol=p)
  }
  
  # Random generation debug
  # rnorm_index <- 1
  # rgig_index <- 1
  # rgamma_index <- 1
  # gamma_vec <- numeric((burnin + nmc) * p)
  # rnorm_vec <- numeric((burnin + nmc) * p * (p - 1))
  # rgig_vec <- numeric((burnin + nmc) * p * (p - 1))
  # runi_vec <- c()
  
  # Begin main MCMC loop
  for (iter in 1:(burnin + nmc)) {
    
    # Gibb's sampler for Omega with Hao-Wang's decomposition
    for (i in 1:p) {
      
      # Selection of relevant indices
      ind_noi <- ind_noi_all[,i]
      
      # Select S column data
      S_21 <- S[ind_noi, i]
      
      # Only S is needed for Wishart
      if (prior != 'Wishart') {
        sigma_11 <- cur_sigma[ind_noi, ind_noi]
        sigma_12 <- cur_sigma[ind_noi, i]
        sigma_22 <- cur_sigma[i, i]
        vec_acc_21 <- matrix_acc_gibbs[ind_noi, i]
        tau_12 <- tau[ind_noi, i]
        inv_omega_11 <- sigma_11 - sigma_12 %*% t(sigma_12) / sigma_22
      }
      else {
        inv_omega_11 <- solve(omega[ind_noi, ind_noi])
      }
      
      if (prior == 'Wishart') {
        shape <- (dof + n - p + 1) / 2
        scale <- 2 / (S[i, i] + 1)
        inv_c <- (S[i, i] + 1) * inv_omega_11
        mu_i <- -solve(inv_c, S_21)
      }
      else if (prior == 'BGL') {
        shape <- n / 2 + 1
        scale <- 2 / (S[i, i] + lambda)
        inv_c <- (
          diag(1 / tau_12, nrow=(p - 1)) + (S[i, i] + lambda) * inv_omega_11
        )
        mu_i <- -solve(
          inv_c, (S_21 + (vec_acc_21 / tau_12))
        )
      }
      else if (prior == 'GHS') {
        nu_12 <- nu[ind_noi, i]
        shape <- n / 2 + 1
        scale <- 2 / (S[i, i] + 1 / lambda)
        inv_c <- (
          diag(1 / (tau_12 * lambda^2), nrow=(p - 1)) + 
          (S[i, i] + 1 / lambda) * inv_omega_11
        )
        mu_i <- -solve(
          inv_c, (S_21 + (vec_acc_21 / (tau_12 * lambda^2)))
        )
      }
      
      # If burnin is complete, save 
      if ((iter > burnin) & (i == p)) {
        inv_c_store[, , (iter - burnin)] <- inv_c
        mean_vec_store[, (iter - burnin)] <- mu_i
      }
      
      # Sampling from the Normal density of Equation (15) in the paper
      cur_rnorm <- rnorm(p - 1)
      # rnorm_vec[rnorm_index:(rnorm_index + p - 2)] <- cur_rnorm
      # rnorm_index <- rnorm_index + p - 1
      beta <- mu_i + solve(chol(inv_c), cur_rnorm)
      
      # Sampling from the gamma density of EQ 15
      gamma_param <- rgamma(1, shape, 1 / scale)
      # gamma_vec[rgamma_index] <- gamma_param
      # rgamma_index <- rgamma_index + 1
      
      omega_22 <- gamma_param + t(beta) %*% inv_omega_11 %*% beta
      
      # Update omega
      omega[i, ind_noi] <- beta
      omega[ind_noi, i] <- beta
      omega[i, i] <- omega_22
      
      # cat("R code omega iter:", iter, "ith: ", i, "\n")
      # print(omega)
      
      if (prior != 'Wishart') {
        # Update sigma
        omega_beta <- inv_omega_11 %*% beta;
        sigma_11 <- inv_omega_11 + (omega_beta %*% t(omega_beta) / gamma_param)
        sigma_12 <- -omega_beta / gamma_param
        sigma_22 <- 1 / gamma_param
        cur_sigma[ind_noi, ind_noi] <- sigma_11
        cur_sigma[ind_noi, i] <- sigma_12
        cur_sigma[i, ind_noi] <- sigma_12
        cur_sigma[i, i] <- sigma_22
        
        # cat("R code cur_sigma iter:", iter, "ith: ", i, "\n")
        # print(cur_sigma)
        
        # Update tau and nu if needed
        if (prior == 'BGL') {
          mu_prime_sq <- lambda^2 / (beta + vec_acc_21)^2
          a_gig_tau <- lambda^2 / mu_prime_sq
          u_12 <- numeric(p - 1)
          for (tau_id in 1:(p - 1)){
            # cat("R code ith:", i, "iter:", tau_id, "\n")
            # print(a_gig_tau[tau_id])
            debug_ls <- gigrnd(-1 / 2, a_gig_tau[tau_id], lambda^2, 1)
            u_12[tau_id] <- debug_ls
            # cat("R code sampled tau: ", 1 / debug_ls$ans, "\n")
            # runi_vec <- c(runi_vec, debug_ls$usample)
            # rgig_vec[rgig_index] <- u_12[tau_id]
            # rgig_index <- rgig_index + 1
          }
          tau_12 <- 1 / u_12
        }
        else if (prior == 'GHS') {
          rate <- beta^2 / (2 * lambda^2) + 1 / nu_12
          tau_12 <- 1 / rgamma(length(rate), 1, rate)
          nu_12 <- 1 / rgamma(length(rate), 1, 1 + 1 / tau_12)
          nu[i, ind_noi] <- nu_12
          nu[ind_noi, i] <- nu_12
        }
        
        tau[i, ind_noi] <- tau_12
        tau[ind_noi, i] <- tau_12
        
        # cat("R code tau iter:", iter, "ith: ", i, "\n")
        # print(tau)
        
      }
    }
    
    if (iter > burnin) {
      omega_acc <- omega_acc + omega
      tau_acc <- tau_acc + tau
    }
  }
  
  ind_noi <- ind_noi_all[, p]
  omega_acc <- omega_acc / nmc
  tau_acc <- tau_acc / nmc
  
  # Calc eq 9
  for (sample_index in 1:nmc) {
    normal_density_acc <- normal_density_acc + mvtnorm::dmvnorm(
      omega_acc[p, ind_noi], mean_vec_store[, sample_index], 
      solve(inv_c_store[, , sample_index])
    )
  }
  MC_avg_eq_9 <- log(normal_density_acc / nmc)
  
  ### R code time profiling ###
  calc_time <- proc.time() - start_time_mcmc_hw
  g_time_env$mcmc_hw_calc_time <- (
    g_time_env$mcmc_hw_calc_time + calc_time
  )
  cat("HAO WANG CALC TIME: \n")
  print(calc_time)
  ######################
  
  return(list(
    MC_avg_eq_9=MC_avg_eq_9,
    post_mean_omega=omega_acc,
    post_mean_tau=tau_acc
  ))
}