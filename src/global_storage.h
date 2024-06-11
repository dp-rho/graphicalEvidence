/* global_storage.h */

/* Global constants for prior coding and size */
#define WISHART (0)
#define BGL (1)
#define GHS (2)
#define MAX_DIM (100)

/* Global info storage/args for manual LAPACK/BLAS calls  */
extern int info_int;
extern char uplo;
extern double one;
extern int nrhs;

/* Global static storage */
extern int g_ipv[MAX_DIM];
extern double g_vec1[MAX_DIM];
extern double g_vec2[MAX_DIM];
extern double g_mat1[MAX_DIM * MAX_DIM];
