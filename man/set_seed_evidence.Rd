% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/set_seed_evidence.R
\name{set_seed_evidence}
\alias{set_seed_evidence}
\title{Set the Random Seed}
\usage{
set_seed_evidence(seed)
}
\arguments{
\item{seed}{a random seed that will be passed to the interpreted random 
number generator using set.seed, and will be passed to the compiled random
number generator using private Rcpp package function set_seed}
}
\description{
Sets the random seed of both the R session (using set.seed) and the compiled
sampler, as both samplers are used during any calls to evidence(...) or 
prior_sampling(...)
}
\examples{
set_seed_evidence(42)
}
