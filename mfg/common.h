#ifndef MFG_COMMON_H
#define MFG_COMMON_H

#include "consts.h"

double func_u(double t, double x);

double *solve_1(int *grid, int *gridPr);
double *calc_error_1(int *grid, double *solution, double tt, int nx_1, double hx);
double *calc_exact_1(int *grid, double t, int nx_1, double hx);

#endif //MFG_COMMON_H