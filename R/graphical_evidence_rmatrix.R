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
  lambda = 0,
  alpha = 0,
  V = 0,
  print_progress = FALSE
) {
  
  # Ensure scale matrix is matrix type
  V <- as.matrix(V)
  
  # Initialize storage for I_{p-j+1} -IV_{p-j+1} for every j from 1 to p
  log_ratio_of_likelihood <- numeric(p - 1)
  
  # Errors possible when calculating log-marginal -> log f(y_{1:p}),
  # return NA in that case
  result <- tryCatch({
    
    # storage of posterior means
    last_col_store <- vector("list", length = p)
    
    # Accumulator for linear shifts
    matrix_acc <- matrix(numeric(p * p), nrow=p)
    
    # Main loop from 1 to p, for Wishart only this can be run in parallel
    if (prior != 'Wishart') { 
      for (num_rmat in 1:(p - 1)) {
        
        updated_density <- graphical_evidence_rmatrix_iteration(
          xx, S, n, p, burnin, nmc, prior, lambda, alpha, V, print_progress,
          matrix_acc, num_rmat
        )
        
        log_ratio_of_likelihood[num_rmat] <- updated_density$log_ratio_of_likelihood
        last_col_store[[num_rmat]] <- updated_density$last_col_store_item
        matrix_acc <- updated_density$matrix_acc
      }

      # Calculate final iteration
      direct_eval_log_prior_density <- graphical_evidence_rmatrix_last_iter(
        xx, S, n, p, nmc, prior, lambda, alpha, V, last_col_store, matrix_acc
      )
      
      return(sum(log_ratio_of_likelihood) + direct_eval_log_prior_density)
    }
    
    # Use doParallel to execute top level loop for Wishart case
    else {
      par_res <- foreach::foreach(num_rmat = 1:(p - 1)) %dopar% {
        
        graphical_evidence_rmatrix_iteration(
          xx, S, n, p, burnin, nmc, prior, lambda, alpha, V, print_progress,
          matrix_acc, num_rmat
        )
      }
      
      # Aggregate parallel computation results
      for (i in 1:(p - 1)) {
        last_col_store[[i]] <- par_res[[i]]$last_col_store_item
        log_ratio_of_likelihood[i] <- par_res[[i]]$log_ratio_of_likelihood
      }
      
      # Calculate final iteration
      direct_eval_log_prior_density <- graphical_evidence_rmatrix_last_iter(
        xx, S, n, p, nmc, prior, lambda, alpha, V, last_col_store, matrix_acc
      )

      return(sum(log_ratio_of_likelihood) + direct_eval_log_prior_density)
    }
    
  }, error = function(e) {
    print(e)
    return(NaN)
  })
  
  return(result)
}