/* graphical_evidence.h */

#pragma once

/* Parallel execution enabled through OpenMP  */
// [[Rcpp:NOT:plugins(openmp)]]
#include <omp.h>

/* RcppArmadillo used for some linear algebra structures  */
// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>

/* Package headers  */
#include <iostream>
#include <chrono>
#include <random>
#include <cmath>
#include <lapacke.h>
#include <cblas.h>
#include "simd_intrinsics.h"
#include "global_storage.h"
#include "FunctionTimer.h"
#include "GammaSampler.h"
#include "prototypes.h"
#include "looping_dmvnrm_arma.h"
#include "inject_random.h"

using namespace Rcpp;