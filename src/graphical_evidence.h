/* graphical_evidence.h */

#pragma once

/* Parallel execution enabled through OpenMP  */
// [[Rcpp::plugins(openmp)]]
#include <omp.h>

/* RcppArmadillo used as wrapper for LAPACK/BLAS routines */
// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>

/* Package headers  */
#include <iostream>
#include <chrono>
#include "FunctionTimer.h"
#include "prototypes.h"
#include "log_dmvnrm_arma_vec.h"

using namespace Rcpp;