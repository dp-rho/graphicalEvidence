#include "graphical_evidence.h"

/* Global random variables injected from R session  */
double* g_rgamma_vec;
double* g_rnorm_vec;
double* g_rgig_vec;
double* g_runi_vec;

/* Count variables to track index for globabl random vectors  */
unsigned int g_rgamma_index;
unsigned int g_rnorm_index;
unsigned int g_rgig_index;
unsigned int g_runi_index;


/*
 * Binds random variables generated in R session to compiled code
 */
// [[Rcpp::export]]
void bind_random_samples(
  NumericVector rgammas,
  NumericVector rnorms
) {

  /* Init indices */
  g_rgamma_index = 0;
  g_rnorm_index = 0;

  /* Bind gamma random variales */
  if (g_rgamma_vec) free(g_rgamma_vec);
  g_rgamma_vec = (double*) malloc(sizeof(double) * rgammas.size());
  memcpy(g_rgamma_vec, rgammas.begin(), sizeof(double) * rgammas.size());

  /* Bind normal random variables */
  if (g_rnorm_vec) free(g_rnorm_vec);
  g_rnorm_vec = (double*) malloc(sizeof(double) * rnorms.size());
  memcpy(g_rnorm_vec, rnorms.begin(), sizeof(double) * rnorms.size());
}

/*
 * Binds random variables generated in R session to compiled code
 */
 // [[Rcpp::export]]
void bind_random_samples_rmatrix(
  NumericVector rgammas,
  NumericVector rnorms,
  NumericVector rgigs,
  NumericVector runis
) {

  /* Init indices */
  g_rgamma_index = 0;
  g_rnorm_index = 0;

  /* Bind gamma random variales */
  if (g_rgamma_vec) free(g_rgamma_vec);
  g_rgamma_vec = (double*)malloc(sizeof(double) * rgammas.size());
  memcpy(g_rgamma_vec, rgammas.begin(), sizeof(double) * rgammas.size());

  /* Bind normal random variables */
  if (g_rnorm_vec) free(g_rnorm_vec);
  g_rnorm_vec = (double*)malloc(sizeof(double) * rnorms.size());
  memcpy(g_rnorm_vec, rnorms.begin(), sizeof(double) * rnorms.size());

  /* Bind generalized inverse gaussian random variables */
  if (g_rgig_vec) free(g_rgig_vec);
  g_rgig_vec = (double*)malloc(sizeof(double) * rgigs.size());
  memcpy(g_rgig_vec, rgigs.begin(), sizeof(double) * rgigs.size());

  /* Bind uniform random variables */
  if (g_runi_vec) free(g_runi_vec);
  g_runi_vec = (double*)malloc(sizeof(double) * runis.size());
  memcpy(g_runi_vec, runis.begin(), sizeof(double) * runis.size());
}


/*
 * Extracts random sampled values from global rnorm values and updates index
 */

void extract_rnorm(double* mem_to_write, unsigned int sample_size) {

  /* Copy random samples and update index */
  memcpy(mem_to_write, g_rnorm_vec + g_rnorm_index, sizeof(double) * sample_size);
  g_rnorm_index += sample_size;
}


/*
 * Extracts random sampled values from global rgamma values and updates index
 */

double extract_rgamma() {
  return g_rgamma_vec[g_rgamma_index++];
}


/*
 * Extracts random sampled values from global rgig values and updates index
 */

double extract_rgig() {
  return g_rgig_vec[g_rgig_index++];
}


/*
 * Extracts random sampled values from global runi values and updates index
 */

double extract_runi() {
  return g_runi_vec[g_runi_index++];
}