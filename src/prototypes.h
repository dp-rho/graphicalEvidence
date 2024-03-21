/* prototypes.h */

/* Functions in inject_random.cpp */
double extract_rgamma();
void extract_rnorm(double*, unsigned int);

/* Functions in GammaSampler.cpp  */
double gamma_density(const double, const double, const double);

/* Functions in double looping_dmvnrm_arma.cpp */
double looping_dmvnrm_arma(
    arma::rowvec const&, arma::mat const&, arma::cube const&
);

void looping_mvpdf_process_iteration(
  arma::vec&, arma::rowvec const&, arma::cube const&, arma::mat const&, 
  arma::mat&, arma::rowvec&, arma::mat&, double const, arma::uword
);

void inplace_tri_mat_mult(
    arma::rowvec&, arma::mat const&
);

/* Functions in calc_eq_9.cpp */
double calc_eq_9(
  arma::uvec const&, arma::cube const&, arma::mat const&,
  arma::mat const&, arma::uword const, const unsigned int, 
  const unsigned int
);

/* Functions in solve_mu_reduced_hw.cpp */
arma::vec solve_mu_reduced_hw(
  const unsigned int, arma::uvec const&, arma::uvec const&,
  arma::mat const&, arma::rowvec&, arma::vec&
);

/* Functions in sample_omega_hw.cpp */
void sample_omega_hw(
  const int, const int, const int, const int, arma::vec&,
  arma::mat&, arma::mat&, arma::mat&, arma::mat&, arma::mat&,
  arma::cube&, arma::mat const&, arma::mat const&, arma::umat const&,
  std::vector<arma::uvec> const&, std::vector<arma::uvec> const&,
  arma::mat const&, arma::mat const&
);

/* Functions in initialize_indices.cpp  */
void initialize_indices(
  arma::mat const&, arma::umat&, std::vector<arma::uvec>&,
  std::vector<arma::uvec>&
);

/* Functions in sample_omega_last_col.cpp */
void sample_omega_last_col(
  const unsigned int, const double, const double*, const double,
  arma::vec&, arma::mat&, arma::mat&, arma::vec&, arma::vec&,
  std::vector<arma::uvec> const&, std::vector<arma::uvec> const&,
  arma::umat const&, arma::mat const&, arma::mat const&, arma::mat const&,
  arma::mat const&
);