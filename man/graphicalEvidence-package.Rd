\name{graphicalEvidence-package}
\alias{graphicalEvidence-package}
\alias{graphicalEvidence}
\docType{package}
\title{
  Compute Marginal Likelihood in Gaussian Graphical Models
}
\description{
  This package allows estimation of marginal likelihood in Gaussian graphical 
  models through a novel telescoping block decomposition of the precision 
  matrix. This package also provides an MCMC prior sampler for the priors of 
  BGL, GHS, and G-Wishart.
}
\details{
  This package implements marginal estimation for four priors, "Wishart"", 
  Bayesian Graphical Lasso ("BGL"), graphical horseshoe ("GHS"), and 
  "G-Wishart". An MCMC prior sampler is also provided for "BGL", "GHS", and
  "G-Wishart".
  
  This version was installed directly from GitHub and may be less portable.
  This package attempts to compile and link with C/C++ libraries of OpenMp, 
  OpenBlas, AVX/SSE, and LAPACK directly. A version which uses only 
  RcppArmadillo is available on CRAN using:
  \code{install.packages('graphicalEvidence')}
}
\author{
Maintainer: David Rowe <david@rowe-stats.com>
}
\references{
  Bhadra, A., Sagar, K., Rowe, D., Banerjee, S., & Datta, J. (2022). Evidence 
  Estimation in Gaussian Graphical Models Using a Telescoping Block 
  Decomposition of the Precision Matrix. arXiv preprint arXiv:2205.01016.
}
\keyword{ package }
\seealso{
  \code{\link{test_evidence}}: For basic example of functionality
  \code{\link{evidence}}: For top level estimation function
  \code{\link{prior_sampling}}: For the prior sampler function
}
\examples{
  \dontrun{
     test_results <- test_evidence(num_runs=10, prior_name="Wishart")   
  }
}
