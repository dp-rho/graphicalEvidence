
# R package documentation needed for each parameter

# Top level function that is called by end user
evidence <- function(
  xx,
  S,
  n,
  p,
  burnin,
  nmc,
  prior_name = c('Wishart', 'BGL', 'GHS', 'G_Wishart'),
  alpha = NULL,
  lambda = NULL,
  V = NULL,
  G = NULL
  # Do we need banded_param for this function?
) {
  
  # Match arg on prior name
  prior_name <- match.arg(prior_name)
  
  result <- switch(
    prior_name,
    
    # Not implemented
    'Wishart' = NULL,
    
    # Not implemented
    'BGL' = NULL,
    
    # Not implemented
    'GHS' = NULL,
    
    # First place to start work
    'G_Wishart' = graphical_evidence_G_Wishart(
      xx, S, n, p, burnin, nmc, alpha, V, G
    )
  )
  
  return(result)
}
