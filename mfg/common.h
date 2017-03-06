#ifndef MFG_COMMON_H
#define MFG_COMMON_H

#include "consts.h"

double *solve_1(int *grid, int *gridPr);
double *calc_exact_1(int *grid, double t, int nx_1, double hx);

#endif //MFG_COMMON_H