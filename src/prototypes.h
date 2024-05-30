/* prototypes.h */

/* Functions in inject_random.cpp */
double extract_rgamma();
double extract_rgig();
double extract_runi();
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

void solve_mu_reduced_hw_in_place(
  const unsigned int, arma::uvec const&, arma::uvec const&,
  arma::mat const&, arma::rowvec const&
);

/* Functions in sample_omega_hw.cpp */
void sample_omega_hw(
  const int, const int, const int, const int, arma::vec&,
  arma::mat&, arma::mat&, arma::mat&, arma::mat&, arma::mat&,
  arma::cube&, arma::mat const&, arma::mat const&, arma::umat const&,
  std::vector<arma::uvec> const&, std::vector<arma::uvec> const&,
  arma::mat const&, arma::mat const&
);

/* Functions in sample_omega_hw_rmatrix */
void sample_omega_hw_rmatrix(
  const int, const int, const int, const int, const int, const double,
  arma::vec&, arma::mat&, arma::mat&, arma::mat&, arma::mat&, arma::mat&,
  arma::mat&, arma::mat&, arma::mat&, arma::mat&, arma::cube&, arma::mat const&,
  arma::umat const&, arma::mat const&, const double, const double*
);

void efficient_inv_omega_11_calc(
  arma::mat&, arma::uvec const&, arma::mat const&,
  const unsigned int, const unsigned int
);

void update_sigma_inplace(
  arma::mat&, arma::mat const&, arma::vec const&, arma::uvec const&,
  const double, const unsigned int, const unsigned int
);

/* Functions in initialize_indices.cpp  */
void initialize_indices(
  arma::umat&
);

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

/* Functions in calc_eq_11.cpp  */
double calc_eq_11(
  const double, const double, const double, const double,
  const unsigned int, const unsigned int, arma::vec const&
);

/* Functions in gigrnd.cpp  */
double gigrnd(
  double, double, double
);

double psi(
  double, double, double
);

double dpsi(
  double, double, double
);

double fun_g(
  double, double, double, double, double
);

/* Functions in mcmc_hw_rmatrix */
/* mcmc_hw_rmatrix() not called internally  */
void get_gamma_params_hw_rmatrix(
  double*, double*, const int, const int, const int,
  const double, arma::mat const&
);