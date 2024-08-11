
graphical_evidence_rmatrix_last_iter <- function(
  xx,
  S,
  n,
  p,
  nmc,
  prior,
  lambda,
  alpha,
  V,
  last_col_store,
  matrix_acc
) {
  
  # Reduce data by one column
  reduced_data_xx <- xx[, 1:1]
  p_reduced <- 1
  S_reduced <- as.matrix(t(reduced_data_xx) %*% reduced_data_xx)
  
  # Prior specific gamma parameters
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
  gamma_samples <- rgamma(
    nmc, shape, 1 / scale
  )
  gamma_mean <- mean(gamma_samples)
  last_col_store_item <- gamma_mean
  
  last_col_store[[p]] <- last_col_store_item
  
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
    log(dnorm(reduced_data_xx, 0, 1 / sqrt(gamma_mean)))
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
        (const / 2) * log(det(as.matrix(diag(1) + S_reduced))) - 
        (const / 2) * log(2) - 
        logmvgamma(const / 2, 1)
    )
  }
  
  # Condense with general case
  log_ratio_of_likelihood <- (
    direct_log_data_density - direct_log_post_density
  )
  
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
    if (p == 2) {
      direct_eval_log_prior_density <- (
        direct_eval_log_prior_density - log(0.67)
      )
    }
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
    if (p == 2) {
      direct_eval_log_prior_density <- (
        direct_eval_log_prior_density - log(0.6446)
      )
    }
  }

  return(direct_eval_log_prior_density + log_ratio_of_likelihood)
}