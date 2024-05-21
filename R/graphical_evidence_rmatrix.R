# Top level function for all priors which use the reconstructed matrix
# method to evaluate log prior density directly, specifically, 
# Wishart, BGL, and GHS
graphical_evidence_rmatrix <- function(
  xx,
  S,
  n,
  p,
  burnin,
  nmc,
  prior,
  lambda = NULL,
  alpha = NULL,
  V = NULL
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
    for (num_rmat in 1:p) {
      
      cat(paste0("Working on ", num_rmat, "th row of the telescoping sum\n"))
      
      # Reduce data by one column
      reduced_data_xx <- xx[, 1:(p - num_rmat + 1)]
      p_reduced <- ncol(reduced_data_xx)
      S_reduced <- as.matrix(t(reduced_data_xx) %*% reduced_data_xx)
      
      # General case
      if (num_rmat < p) {
        
        # Matrix accumulator is not needed in Wishart, will never be used
        matrix_acc_gibbs <- matrix_acc[1:p_reduced, 1:p_reduced]
        
        # Run unrestricted Hao Wang sampler which will
        # be used to approximate the Normal density in the
        # evaluation of the term IV_{p-j+1} - Equation (9)
        hw_results <- rmatrix_Hao_Wang(
          S_reduced, n, burnin, nmc, prior, dof=(alpha + 1 - num_rmat),
          matrix_acc_gibbs=matrix_acc_gibbs, lambda=lambda
        )

        fixed_last_col <- (
          hw_results$post_mean_omega[1:(p_reduced - 1), p_reduced]
        )
        
        last_col_results <- rmatrix_last_col_fixed(
          S_reduced, n, burnin, nmc, fixed_last_col, prior,
          dof=(alpha + 1 - num_rmat), matrix_acc_gibbs=matrix_acc_gibbs,
          post_mean_omega=hw_results$post_mean_omega,
          post_mean_tau=hw_results$post_mean_tau, lambda=lambda
        )
        
        # Computing IV_{p-j+1}
        log_normal_posterior_density[num_rmat] <- hw_results$MC_avg_eq_9
        log_gamma_posterior_density[num_rmat] <- last_col_results$MC_avg_eq_11
        log_posterior_density[num_rmat] <- (
          hw_results$MC_avg_eq_9 + last_col_results$MC_avg_eq_11
        )
        
        # Computing I_{p-j+1}
        st_dev <- sqrt(1 / last_col_results$post_mean_omega_22)
        mean_vec <- -(st_dev^2) * (
          reduced_data_xx[, 1:(p_reduced - 1)] %*% 
          as.matrix(hw_results$post_mean_omega[p_reduced, 1:(p_reduced - 1)])
        )
        log_data_likelihood[num_rmat] <- sum(log(
          dnorm(reduced_data_xx[, p_reduced], mean_vec, st_dev)
        ))
        
        # Computing III_{p-j+1}
        log_ratio_of_likelihood[num_rmat] <- (
          log_data_likelihood[num_rmat] - log_posterior_density[num_rmat]
        )
        
        # Not needed for Wishart
        if (prior != 'Wishart') {
          matrix_acc[1:(p - num_rmat), 1:(p - num_rmat)] <- (
            matrix_acc[1:(p - num_rmat), 1:(p - num_rmat)] +
            (1 / last_col_results$post_mean_omega_22) * (
              fixed_last_col %*% t(fixed_last_col)
            )
          )
        }
        
        last_col_store[[num_rmat]] <- c(
          fixed_last_col, last_col_results$post_mean_omega_22
        )
      }
      
      # Final iteration
      else {
        
        if (prior == 'Wishart') {
          shape <- 0.5 * (alpha + (1 - p) + n)
          scale <- 2 / (1 + S_reduced[1])
        }
        else if (prior == 'BGL') {
          shape <- n / 2 + 1
          scale <- 2 / (lambda + S_reduced[1])
        }
        else if (prior == 'GHS') {
          shape <- n / 2 + 1
          scale <- 2 / (1 / lambda + S_reduced[1])
        }
        
        # To evaluate log f(y_1) or the last row in the telescopic sum. Here
        # we compute I_1 and IV_1 only. II_1=0 because no further columns of the
        # data are left and III_1 is taken care of the joint prior evaluation
        gamma_samples <- rgamma_compiled(
          nmc, shape, 1 / scale
        )
        gamma_mean <- mean(gamma_samples)
        
        last_col_store[[num_rmat]] <- gamma_mean
        
        reconstructed_matrix <- matrix(0, nrow=p, ncol=p)
        for (col_id in 1:p) {
          reconstructed_matrix[1:col_id, col_id] <- last_col_store[[
            p - col_id + 1
          ]]
          reconstructed_matrix[col_id, 1:col_id] <- last_col_store[[
            p - col_id + 1
          ]]
        }
        
        for (col_id in 1:(p - 1)) {
          col_len <- col_id + 1
          col_vec <- last_col_store[[p - col_id]][1:(col_len - 1)]
          diag_ele <- last_col_store[[p - col_id]][col_len]
          reconstructed_matrix[1:col_id, 1:col_id] <- (
            reconstructed_matrix[1:col_id, 1:col_id] +
            (1 / diag_ele) * (col_vec %*% t(col_vec))
          )
        }
        
        direct_log_data_density <- sum(
          log(dnorm(reduced_data_xx, 0, 1 / sqrt(gamma_mean))) # inv
        )
        
        # Log posterior density dependent on prior
        if (prior != 'Wishart') {
          direct_log_post_density <- log(
            dgamma(gamma_mean, n / 2  + 1, scale=(2 / (lambda + S_reduced)))
          )
        }
        else {
          const <- alpha + (1 - p) + n
          direct_log_post_density <- (
            ((const - 2) / 2) * log(gamma_mean) - 
            (S_reduced[1] + 1) * gamma_mean / 2 + 
            (const / 2) * log(det(diag(1) + S_reduced)) - 
            (const / 2) * log(2) - 
            logmvgamma(const / 2, 1)
          )
        }
        
        # Condense with general case
        log_data_likelihood[num_rmat] <- direct_log_data_density
        log_posterior_density[num_rmat] <- direct_log_post_density
        log_ratio_of_likelihood[num_rmat] <- direct_log_data_density - direct_log_post_density
        
        # Direct eval of log prior density
        abs_lower_tri_est_mat <- abs(
          reconstructed_matrix * lower.tri(reconstructed_matrix)
        )
        
        # Direct evaluation is dependent on prior
        if (prior == 'Wishart') {
          direct_eval_log_prior_density <- (
            0.5 * (alpha - p - 1) * log(det(reconstructed_matrix)) - 
            0.5 * sum(diag(reconstructed_matrix)) - 0.5 * alpha * p * log(2) - 
            logmvgamma(alpha / 2, p)
          )
          direct_eval_log_prior_density <- (
            direct_eval_log_prior_density + (n / 2) * log(det(V))
          )
        }
        else if (prior == 'BGL') {
          direct_eval_log_prior_density <- (
            0.5 * (p * (p - 1)) * log(lambda / 2) - lambda * 
            sum(abs_lower_tri_est_mat) + p * log(lambda / 2) - 
            (lambda / 2) * sum(diag(reconstructed_matrix))
          )
        }
        else if (prior == 'GHS') {
          lower_tri_iterators <- 1:(p * (p - 1) / 2)
          log_dawson_vals <- numeric(length(lower_tri_iterators))
          abs_lower_tri_vec <- abs_lower_tri_est_mat[
            which(abs_lower_tri_est_mat != 0)
          ]
          
          # Horseshoe density can also be expressed as a
          # Laplace mixture, this is more stable than drawing
          # samples from a Cauchy distribution
          for (lower_tri_iter in lower_tri_iterators) {
            rate <- abs_lower_tri_vec[lower_tri_iter] / lambda
            rand_sample <- rexp(1e4, rate) / sqrt(2)
            rand_sample_p2 <- rand_sample ^ 2
            rand_sample_p4 <- rand_sample ^ 4
            rand_sample_p6 <- rand_sample ^ 6
            
            dawson_num <- (
                1 + (33 / 232) * rand_sample_p2 + (19 / 632) * rand_sample_p4 + 
                (23 / 1471) * rand_sample_p6
            )
            
            dawson_denom <- (
              1 + (517 / 646) * rand_sample_p2 + 
              (58 / 173) * rand_sample_p4 +
              (11 / 262) * rand_sample_p6 +
              (46 / 1471) * (rand_sample ^ 8)
            )
            
            log_dawson_vals[lower_tri_iter] <- log(
              mean(rand_sample * dawson_num / dawson_denom)
            )
          }
          
          direct_eval_log_prior_density <- (
            p * (p - 1) / 2 * (log(2) - 1.5 * log(pi)) + 
            sum(log_dawson_vals) - sum(log(abs_lower_tri_vec)) + 
            p * log(1 / (2 * lambda)) - 
            sum(diag(reconstructed_matrix)) / (2 * lambda)
          )
        }
      }
      
    }
    return(sum(log_ratio_of_likelihood) + direct_eval_log_prior_density)
    
  }, error = function(e) {
    print(e)
    return(NaN)
  })
  
  return(result)
}