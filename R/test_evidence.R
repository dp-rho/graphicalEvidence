#' @title Test Graphical Evidence
#' 
#' @description
#' Tests any of the allowed priors with preexisting test inputs which
#' should yield reproducible results, as the random seed is reset during
#' this function call
#' 
#' @param num_runs An integer number that specifies how many runs of 
#' graphical evidence will be executed on the test parameters, as multiple
#' runs allows us to quantify uncertainty on the estimator.
#' 
#' @param prior_name The name of the prior for being tested with preexisting
#' test parameters, this is one of 'Wishart', 'BGL', 'GHS', 'G_Wishart'
#' 
#' @returns A list of results which contains the mean marginal likelihood, the
#' standard deviation of the estimator, and the raw results in a vector
#' 
#' @examples
#' # Compute the marginal 10 times with random column permutations of the 
#' # preexisting test parameters for G-Wishart prior 
#' test_evidence(10, 'G_Wishart')
#' @export
test_evidence <- function(
    num_runs, 
    prior_name = c('Wishart', 'BGL', 'GHS', 'G_Wishart')    
) {
  
  # Match arg on prior name
  prior_name <- match.arg(prior_name)
  
  # Get predetermined test inputs
  params <- gen_params_evidence(prior_name)
  
  cat("Params are: \n")
  print(params)
  
  # Make results repeatable
  set_seed_evidence(123456789)
  
  # Save start in R program time
  our_time <- proc.time()

  # Call top level function
  estimated_marginal_store <- switch(
    prior_name,
    
    'Wishart' = evidence(
      xx=params$x_mat, burnin=1e3, nmc=5e3, prior_name=prior_name, 
      runs=num_runs, alpha=7, V=params$scale_mat,
    ),
    
    'BGL' = evidence(
      xx=params$x_mat, burnin=1e3, nmc=5e3, prior_name=prior_name, 
      runs=num_runs, lambda=1
    ),
    
    'GHS' = evidence(
      xx=params$x_mat, burnin=1e3, nmc=5e3, prior_name=prior_name, 
      runs=num_runs, lambda=1
    ),
    
    'G_Wishart' = evidence(
      xx=params$x_mat, burnin=2e3, nmc=1e4, prior_name=prior_name, 
      runs=num_runs, alpha=2, V=params$scale_mat, G=params$g_mat
    )
  )
  
  # Cumulative execution time in R program
  cat('Execution time in R program per run (seconds):\n')
  print((proc.time() - our_time) / num_runs)
  
  # Visually plot multiple runs in a histogram
  if (num_runs > 1) {
    hist(
      estimated_marginal_store$results, xlab='estimated marginal', 
      main='Histogram of Results', breaks=40
    )
  }
  return(estimated_marginal_store)
}