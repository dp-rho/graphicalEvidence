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
#include <random>
#include <cmath>
#include <lapacke.h>
#include <cblas.h>
#include <immintrin.h>
#include "global_storage.h"
#include "FunctionTimer.h"
#include "GammaSampler.h"
#include "prototypes.h"
#include "looping_dmvnrm_arma.h"
#include "inject_random.h"

using namespace Rcpp;