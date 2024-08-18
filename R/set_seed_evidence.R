#' set_seed_evidence
#'
#' This is a private function that sets the compiled and interpreted seed
#' so that estimation and sampling done with graphicalEvidence will be 
#' reproducible
#' 
#' @param seed a random seed that will be passed to the interpreted random 
#' number generator using set.seed, and will be passed to the compiled random
#' number generator using private Rcpp package function set_seed
#' 
#' @returns NULL
#' @keywords internal
#' @noRd
set_seed_evidence <- function(
    seed
) {
  
  # Set interpreted seed
  set.seed(seed)
  
  # Set compiled seed
  set_seed(seed)
}