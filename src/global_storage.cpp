#include "graphical_evidence.h"

/* Global info storage for manual LAPACK calls  */
int info_int;

/* Define global static storage */
int g_ipv[MAX_DIM];
double g_vec1[MAX_DIM];
double g_vec2[MAX_DIM];
double g_mat1[MAX_DIM * MAX_DIM];