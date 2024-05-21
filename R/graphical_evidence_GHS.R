
graphical_evidence_GHS <- function(
  xx,
  S,
  n,
  p,
  burnin,
  nmc,
  lambda   
) {
  
  # Initialize storage for I_{p-j+1} -IV_{p-j+1} for every j from 1 to p
  log_ratio_of_likelihood <- numeric(p)
  log_data_likelihood <- numeric(p)
  log_posterior_density <- numeric(p)
  log_normal_posterior_density <- numeric(p)
  log_gamma_posterior_density <- numeric(p)
  
  # Errors possible when calculating log-marginal -> log f(y_{1:p}),
  # return NA in that case
  result <- tryCatch({
    
    # storage of posterior means
    last_col_store <- list()
    
    # Accumulator for linear shifts
    matrix_acc <- matrix(numeric(p * p), nrow=p)
    
    # Main BGL loop from 1 to p
    for (num_BGL in 1:p) {
      
      cat(paste0("Working on ", num_BGL, "th row of the telescoping sum\n"))
      
      # Reduce data by one column
      reduced_data_xx <- xx[, 1:(p - num_BGL + 1)]
      p_reduced <- ncol(reduced_data_xx)
      S_reduced <- as.matrix(t(reduced_data_xx) %*% reduced_data_xx)
      
      # General case
      if (num_BGL < p) {
        
        # Reduce data
        matrix_acc_gibbs <- matrix_acc[1:p_reduced, 1:p_reduced]
        
        # Run unrestricted Hao Wang sampler which will
        # be used to approximate the Normal density in the
        # evaluation of the term IV_{p-j+1} - Equation (9)
        hw_results <- GHS_Hao_Wang(
          S_reduced, n, burnin, nmc, lambda, matrix_acc_gibbs
        )
        
        fixed_last_col <- (
          hw_results$post_mean_omega[1:(p_reduced - 1), p_reduced]
        )
        
        last_col_results <- GHS_last_col_fixed(
          S_reduced, n, burnin, nmc, lambda, fixed_last_col,
          matrix_acc_gibbs, hw_results$post_mean_omega,
          hw_results$post_mean_tau
        )
        
        # Computing IV_{p-j+1}
        log_normal_posterior_density[num_BGL] <- hw_results$MC_avg_eq_9
        log_gamma_posterior_density[num_BGL] <- last_col_results$MC_avg_eq_11
        log_posterior_density[num_BGL] <- (
          hw_results$MC_avg_eq_9 + last_col_results$MC_avg_eq_11
        )
        
        # Computing I_{p-j+1}
        st_dev <- sqrt(1 / last_col_results$post_mean_omega_22)
        mean_vec <- -(st_dev^2) * (
          reduced_data_xx[, 1:(p_reduced - 1)] %*% 
          as.matrix(hw_results$post_mean_omega[p_reduced, 1:(p_reduced - 1)])
        )
        log_data_likelihood[num_BGL] <- sum(log(
          dnorm(reduced_data_xx[, p_reduced], mean_vec, st_dev)
        ))
        
        # Computing III_{p-j+1}
        log_ratio_of_likelihood[num_BGL] <- (
          log_data_likelihood[num_BGL] - log_posterior_density[num_BGL]
        )
        
        matrix_acc[1:(p - num_BGL), 1:(p - num_BGL)] <- (
          matrix_acc[1:(p - num_BGL), 1:(p - num_BGL)] +
          (1 / last_col_results$post_mean_omega_22) * (
            fixed_last_col %*% t(fixed_last_col)
          )
        )
        last_col_store[[num_BGL]] <- c(
          fixed_last_col, last_col_results$post_mean_omega_22
        )
      }
      
      # Final iteration
      else {
        
        # To evaluate log f(y_1) or the last row in the telescopic sum. Here
        # we compute I_1 and IV_1 only. II_1=0 because no further columns of the
        # data are left and III_1 is taken care of the joint prior evaluation
        gamma_samples <- rgamma(
          nmc, n / 2 + 1, scale=(2 / (lambda + S_reduced[1]))
        )
        gamma_mean <- mean(gamma_samples)
        
        last_col_store[[num_BGL]] <- gamma_mean
        
        reconstructed_matrix <- matrix(0, nrow=p, ncol=p)
        for (col_id in 1:p) {
          # browser()
          reconstructed_matrix[1:col_id, col_id] <- last_col_store[[
            p - col_id + 1
          ]]
          reconstructed_matrix[col_id, 1:col_id] <- last_col_store[[
            p - col_id + 1
          ]]
        }
        
        for (col_id in 1:(p - 1)) {
          # browser()
          col_len <- col_id + 1
          col_vec <- last_col_store[[p - col_id]][1:(col_len - 1)]
          diag_ele <- last_col_store[[p - col_id]][col_len]
          reconstructed_matrix[1:col_id, 1:col_id] <- (
            reconstructed_matrix[1:col_id, 1:col_id] +
            (1 / diag_ele) * (col_vec %*% t(col_vec))
          )
        }
        
        direct_log_data_density <- sum(log(dnorm(reduced_data_xx, 0, solve(sqrt(gamma_mean)))))
        direct_log_post_density <- log(dgamma(gamma_mean, n / 2  + 1, scale=(2 / (lambda + S_reduced))))
        
        # Condense with general case
        log_data_likelihood[num_BGL] <- direct_log_data_density
        log_posterior_density[num_BGL] <- direct_log_post_density
        log_ratio_of_likelihood[num_BGL] <- direct_log_data_density - direct_log_post_density
        
        direct_eval_log_prior_density <- (
          0.5 * (p * (p - 1)) * log(lambda / 2) - lambda * 
          sum(abs(
            reconstructed_matrix * lower.tri(reconstructed_matrix)
          )) + p * log(lambda / 2) - (lambda / 2) * sum(diag(reconstructed_matrix))
        )
        # browser()
          
      }
    }
    return(sum(log_ratio_of_likelihood) + direct_eval_log_prior_density)
    
  }, error = function(e) {
    # print(e)
    return(NaN)
  })
  print(result)
  return(result)
}