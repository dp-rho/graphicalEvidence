/* global_storage.h */

#define MAX_DIM (100)

/* Global info storage for manual LAPACK calls  */
extern int info_int;

/* Global static storage */
extern int g_ipv[MAX_DIM];
extern double g_vec1[MAX_DIM];
extern double g_vec2[MAX_DIM];
extern double g_mat1[MAX_DIM * MAX_DIM];
