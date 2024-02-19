/* prototypes.h */

/* Functions in log_dmvnrm_arma_vec.cpp */
arma::vec log_dmvnrm_arma_vec(
    arma::rowvec const&, arma::mat const&, arma::cube const&
);

void looping_mvpdf_process_iteration(
  arma::vec&, arma::rowvec const&, arma::cube const&, arma::mat const&, 
  arma::rowvec&, arma::mat&, double const, arma::uword
);

void inplace_tri_mat_mult(
    arma::rowvec&, arma::mat const&
);

/* Functions in calc_eq_9.cpp */
Rcpp::List calc_eq_9(
    Rcpp::NumericVector, Rcpp::NumericVector, Rcpp::NumericVector,
    Rcpp::NumericVector, int, int
);