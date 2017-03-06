#ifndef MFG_COMMON_H
#define MFG_COMMON_H

#include "consts.h"

double *get_exact_solution(double *vec, int nx_1, double time, double hx, double a,
                           double (*analytical_solution)(double, double));

double *solve_1(double *vec, double *vecPr);
double analytical_solution_1(double time, double dot);
double *calc_exact_1(int *vec, int nx_1, double time, double hx);

#endif //MFG_COMMON_H