
# R package documentation needed for each parameter

# Top level function that is called by end user
evidence <- function(
  xx,
  burnin,
  nmc,
  prior_name = c('Wishart', 'BGL', 'GHS', 'G_Wishart'),
  runs = 1,
  print_progress = FALSE,
  alpha = NULL,
  lambda = NULL,
  V = NULL,
  G = NULL
) {
  
  # Match arg on prior name
  prior_name <- match.arg(prior_name)
  
  # Detect available cores and assign to compiled program
  num_cores <- parallel::detectCores()
  set_cores(num_cores)
  
  # Extract sample size and dimension of xx
  xx <- as.matrix(xx)
  p <- ncol(xx)
  n <- nrow(xx)
  
  # Generate permutations for columns of xx
  permutation_matrix <- t(replicate(runs, sample(p)))
  
  # Store repeated results 
  results <- numeric(runs)
  
  for (i in 1:runs) {
    
    # Use a random permutation of columns of xx
    test_perm <- permutation_matrix[i, ]
    
    # Permute xx columns and then calculate S for multiple runs
    if (i == 1) {
      
      xx_perm <- xx
      if (!is.null(G)) {
        G_perm <- G 
      }
      if (!is.null(V)) {
        V_perm <- V 
      }
    }
    else {
      
      xx_perm <- as.matrix(xx[, test_perm])
      if (!is.null(G)) {
        G_perm <- as.matrix(G[test_perm, test_perm])
      }
      if (!is.null(V)) {
        V_perm <- as.matrix(V[test_perm, test_perm])
      }
    }
    
    # Calculate sample covariance matrix
    S <- as.matrix(t(xx_perm) %*% xx_perm)

    # Run graphical evidence method on identified prior
    results[i] <- switch(
      prior_name,
      
      # Initial R implementation
      'Wishart' = graphical_evidence_rmatrix(
        xx_perm, S, n, p, burnin, nmc, prior_name, alpha=alpha, V=V_perm,
        print_progress=print_progress
      ),
      
      # Largely implemented in C++
      'BGL' = graphical_evidence_rmatrix(
        xx_perm, S, n, p, burnin, nmc, prior_name, lambda=lambda,
        print_progress=print_progress
      ),
      
      # Largely implemented in C++
      'GHS' = graphical_evidence_rmatrix(
        xx_perm, S, n, p, burnin, nmc, prior_name, lambda=lambda,
        print_progress=print_progress
      ),
      
      # Largely implemented in C++
      'G_Wishart' = graphical_evidence_G_Wishart(
        xx_perm, S, n, p, burnin, nmc, alpha, V_perm, G_perm,
        print_progress=print_progress
      )
    )
    if (print_progress) {
      cat("done with run", i, '\n')
    }
  }
  
  # Filter out infinite and NA values from results
  init_len <- length(results)
  results <- results[!is.na(results) & !is.infinite(results)]
  trunc_len <- length(results)
  cat(init_len - trunc_len, "runs excluded for NaN\n")
  
  return(
    list(mean=mean(results), sd=sqrt(var(results)), results=results)
  )
}
