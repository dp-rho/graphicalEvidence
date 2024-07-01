#include "graphical_evidence.h"

/*
 * Computationally efficient calculation of inverse of omega_11:
 * For this implementation, inverse of omega_11 can be calculated as
 * sigma_11 - sigma_12 %*% t(sigma_12) / sigma_22
 */

void efficient_inv_omega_11_calc_no_simd(
  arma::mat& inv_omega_11,
  arma::uvec const& ind_noi,
  arma::mat const& sigma,
  const unsigned int p,
  const unsigned int ith
) {

  /* Iterate over cols  */
  for (unsigned int i = 0; i < (p - 1); i++) {

    /* Iterate over rows  */
    for (unsigned int j = 0; j < (p - 1); j++) {

      inv_omega_11.at(j, i) = sigma.at(ind_noi[j], ind_noi[i]) - (
        sigma.at(ind_noi[j], ith) * sigma.at(ind_noi[i], ith) /
        sigma.at(ith, ith)
      );
    }
  }
}


/*
 * Computationally efficient calculation of inverse of omega_11:
 * For this implementation, inverse of omega_11 can be calculated as
 * sigma_11 - sigma_12 %*% t(sigma_12) / sigma_22
 * Uses explicit SIMD with support for SSE, AVX, AVX512
 */

void efficient_inv_omega_11_calc(
  arma::mat& inv_omega_11,
  arma::uvec const& ind_noi,
  arma::mat const& sigma,
  const unsigned int p,
  const unsigned int ith
) {

  /* Check for supported SIMD */

#if defined(__AVX512F__) || defined(__AVX__) || defined(__SSE2__)

  /* Given the SIMD capabilities of the CPU, identify the largest group */
  /* of columns and rows that can computed using SIMD before and after  */
  /* the separation at index ith                                        */
  unsigned int u1 = (ith / SIMD_WIDTH) * SIMD_WIDTH;
  unsigned int r1 = (ith % SIMD_WIDTH) + u1;
  unsigned int u2 = ((p - ith - 1) / SIMD_WIDTH) * SIMD_WIDTH + ith + 1;
  unsigned int r2 = ((p - ith - 1) % SIMD_WIDTH) + u2;

  /* The ith column pointer of sigma */
  const double* sigma_col_i = sigma.colptr(ith);

  /* Load sigma[ith, ith] to SIMD register  */
  _simd_type sigma_ii = _simd_set1_pd(sigma.at(ith, ith));

  /* Iterate over cols from [0, ith): */
  for (unsigned int cur_col = 0; cur_col < ith; cur_col++) {

    /* The current column pointers of sigma and inv_omega_11 */
    const double* sigma_cur_col = sigma.colptr(cur_col);
    double* inv11_cur_col = inv_omega_11.colptr(cur_col);

    /* Load sigma[cur_col, ith] to SIMD register  */
    _simd_type m1 = _simd_set1_pd(sigma.at(cur_col, ith));

    /* Iterate over rows from [0, ith): */
    /* First SIMD iteration up to u1    */
    for (unsigned int cur_row = 0; cur_row < u1; cur_row += SIMD_WIDTH) {

      /* Load sigma[cur_row, ith] */
      _simd_type m2 = _simd_loadu_pd(&sigma_col_i[cur_row]);

      /* Multiply sigma[cur_col, ith] * sigma[cur_row, ith] */
      m2 = _simd_mul_pd(m1, m2);

      /* Divide product by sigma[ith, ith]  */
      m2 = _simd_div_pd(m2, sigma_ii);

      /* Load in sigma[cur_row, cur_col]  */
      _simd_type res = _simd_loadu_pd(&sigma_cur_col[cur_row]);

      /* Calculate final result to update in inv_omega_11, note     */
      /* that storing this value is done at cur_row and cur_col     */
      /* relative to sigma for this block, however, this varies     */
      /* for other looping blocks, as inv_omega_11 has dim (p - 1)  */
      res = _simd_sub_pd(res, m2);
      _simd_storeu_pd(&inv11_cur_col[cur_row], res);
    }

    /* Finish remaining rows using standard looping from u1 to r1 */
    for (unsigned int cur_row = u1; cur_row < r1; cur_row++) {
      inv_omega_11.at(cur_row, cur_col) = sigma.at(cur_row, cur_col) - (
        sigma.at(cur_row, ith) * sigma.at(cur_col, ith) /
        sigma.at(ith, ith)
      );
    }

    /* Iterate over rows from [ith + 1, p): */
    /* First SIMD iteration up to u2    */
    for (unsigned int cur_row = ith + 1; cur_row < u2; cur_row += SIMD_WIDTH) {

      /* Load sigma[cur_row, ith] */
      _simd_type m2 = _simd_loadu_pd(&sigma_col_i[cur_row]);

      /* Multiply sigma[cur_col, ith] * sigma[cur_row, ith] */
      m2 = _simd_mul_pd(m1, m2);

      /* Divide product by sigma[ith, ith]  */
      m2 = _simd_div_pd(m2, sigma_ii);

      /* Load in sigma[cur_row, cur_col]  */
      _simd_type res = _simd_loadu_pd(&sigma_cur_col[cur_row]);

      /* Calculate final result to update in inv_omega_11, note     */
      /* that storing this value is now at cur_row - 1 and cur_col  */
      res = _simd_sub_pd(res, m2);
      _simd_storeu_pd(&inv11_cur_col[cur_row - 1], res);
    }

    /* Finish remaining rows using standard looping from u2 to r2 */
    for (unsigned int cur_row = u2; cur_row < r2; cur_row++) {
      inv_omega_11.at(cur_row - 1, cur_col) = sigma.at(cur_row, cur_col) - (
        sigma.at(cur_row, ith) * sigma.at(cur_col, ith) /
        sigma.at(ith, ith)
      );
    }
  }

  /* Iterate over cols from [ith + 1, p)  */
  for (unsigned int cur_col = ith + 1; cur_col < p; cur_col++) {

    /* The current column pointers of sigma and inv_omega_11 */
    const double* sigma_cur_col = sigma.colptr(cur_col);
    double* inv11_cur_col = inv_omega_11.colptr(cur_col - 1);

    /* Load sigma[cur_col, ith] to SIMD register  */
    _simd_type m1 = _simd_set1_pd(sigma.at(cur_col, ith));

    /* Iterate over rows from [0, ith): */
    /* First SIMD iteration up to u1    */
    for (unsigned int cur_row = 0; cur_row < u1; cur_row += SIMD_WIDTH) {

      /* Load sigma[cur_row, ith] */
      _simd_type m2 = _simd_loadu_pd(&sigma_col_i[cur_row]);

      /* Multiply sigma[cur_col, ith] * sigma[cur_row, ith] */
      m2 = _simd_mul_pd(m1, m2);

      /* Divide product by sigma[ith, ith]  */
      m2 = _simd_div_pd(m2, sigma_ii);

      /* Load in sigma[cur_row, cur_col]  */
      _simd_type res = _simd_loadu_pd(&sigma_cur_col[cur_row]);

      /* Calculate final result to update in inv_omega_11, note     */
      /* that storing this value is done at cur_row and cur_col - 1 */
      res = _simd_sub_pd(res, m2);
      _simd_storeu_pd(&inv11_cur_col[cur_row], res);
    }

    /* Finish remaining rows using standard looping from u1 to r1 */
    for (unsigned int cur_row = u1; cur_row < r1; cur_row++) {
      inv_omega_11.at(cur_row, cur_col - 1) = sigma.at(cur_row, cur_col) - (
        sigma.at(cur_row, ith) * sigma.at(cur_col, ith) /
        sigma.at(ith, ith)
      );
    }
    /*
    if (debug) {
      arma::cout << "finished normal row iteration b3\n";
      arma::cout << "cur inv_omega_11: \n" << inv_omega_11 << arma::endl;
    } */

    /* Iterate over rows from [ith + 1, p): */
    /* First SIMD iteration up to u2        */
    for (unsigned int cur_row = ith + 1; cur_row < u2; cur_row += SIMD_WIDTH) {

      /* Load sigma[cur_row, ith] */
      _simd_type m2 = _simd_loadu_pd(&sigma_col_i[cur_row]);

      /* Multiply sigma[cur_col, ith] * sigma[cur_row, ith] */
      m2 = _simd_mul_pd(m1, m2);

      /* Divide product by sigma[ith, ith]  */
      m2 = _simd_div_pd(m2, sigma_ii);

      /* Load in sigma[cur_row, cur_col]  */
      _simd_type res = _simd_loadu_pd(&sigma_cur_col[cur_row]);

      /* Calculate final result to update in inv_omega_11, note         */
      /* that storing this value is now at cur_row - 1 and cur_col - 1  */
      res = _simd_sub_pd(res, m2);
      _simd_storeu_pd(&inv11_cur_col[cur_row - 1], res);
    }

    /* Finish remaining rows using standard looping from u2 to r2 */
    for (unsigned int cur_row = u2; cur_row < r2; cur_row++) {
      inv_omega_11.at(cur_row - 1, cur_col - 1) = sigma.at(cur_row, cur_col) - (
        sigma.at(cur_row, ith) * sigma.at(cur_col, ith) /
        sigma.at(ith, ith)
      );
    }
  }

#else

  /* Same calculation with no explicit SIMD, i.e., one iteration at a time  */
  efficient_inv_omega_11_calc_no_simd(
    inv_omega_11, ind_noi, sigma, p, ith
  );

#endif

}