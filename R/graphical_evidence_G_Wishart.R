
# R wrapper for compiled code of G_Wishart prior implementation,
# eventually can be replaced entirely by a call to Rcpp function
graphical_evidence_G_Wishart(
  xx,
  S,
  n,
  p,
  burnin,
  nmc,
  alpha,
  G,
  banded_param
) {
  
  # Check if banded matrix G is requested
  is_banded <- 0
  if (!is.null(banded_param)) {
    is_banded <- 1
  }
  
  ## Do we actually need is_banded ?  It appears this is just used to 
  ## validate results as there exists a closed form solution when
  ## a banded matrix is passed.  I am thinking we don't need the 
  ## banded_param unless we want to add in the validation to this
  ## function call as well
  
  
  # Call to Rcpp function should be here, might cut out banded inputs:
  # compiled_G_Wishart(
  #   xx, S, n, p, burning, nmc, alpha, G, is_banded, banded_param
  # )
  
  
  
}