#ifndef FEM_CIRCLE_COMMON_H
#define FEM_CIRCLE_COMMON_H

#include "consts.h"

double func_u(double t, double x, double y);
double func_v(double t, double x, double y);

double *solve_1(int *grid, int *gridPr);

double *calc_error_1(int *grid, double *solution, double tt, int nx3_1, int ny3_1, double hx_smallest, double hy_smallest, int max_lvl);

double *calc_exact_1(int *grid, double t, int nx3_1, int ny3_1, double hx_smallest, double hy_smallest, int max_lvl);

#endif //FEM_CIRCLE_COMMON_H