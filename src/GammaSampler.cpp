#include "graphical_evidence.h"

/* Initialize the global gamma sampler  */
GammaSampler g_rgamma;


/*
 * Returns the density of the gamma distribution under the specified 
 * shape and scale parameters for input x
 */

double gamma_density(const double x, const double shape, const double scale) {
  return pow(x, shape - 1) * exp(-x / scale) / (pow(scale, shape) * tgamma(shape));
}