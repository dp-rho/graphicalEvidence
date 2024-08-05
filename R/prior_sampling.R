

# Samples nmc_prior iterations from the specified density after burnin_prior
# initial iterations that are not recorded
prior_sampling <- function(
    p, 
    burnin_prior, 
    nmc_prior,
    prior_name = c('BGL', 'GHS', 'G_Wishart'),
    seed = NULL,
    G_mat_adj = NULL, 
    scale_matrix = NULL, 
    alpha = NULL,
    lambda = NULL
) {
  
  # Match arg on prior name
  prior_name <- match.arg(prior_name)
  
  # Set compiled seed and interpreted seed if input seed is not null
  if (!is.null(seed)) {
    set.seed(seed)
    set_seed(seed)
  }
  
  coded_prior <- switch(
    prior_name,
    'BGL' = 1,
    'GHS'= 2
  )

  # Draw samples using compiled sampler dependent on prior
  sampled_omegas <- switch(
    prior_name,
    
    # Implemented in C++
    'BGL' = prior_sampler_rmatrix(
      p, burnin_prior, nmc_prior, coded_prior, lambda
    ),
    
    # Implemented in C++
    'GHS' = prior_sampler_rmatrix(
      p, burnin_prior, nmc_prior, coded_prior, lambda
    ),
    
    # Implemented in C++
    'G_Wishart' = prior_sampler_G_Wishart(
      p, burnin_prior, nmc_prior, as.matrix(G_mat_adj), 
      as.matrix(scale_matrix), alpha
    )
  )
  
  return(sampled_omegas)
}