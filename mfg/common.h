#ifndef MFG_COMMON_H
#define MFG_COMMON_H

#include "consts.h"

double *get_exact_solution(double *vec, int nx_1, double time, double hx, double a,
                           double (*analytical_solution)(double, double));

double *solve_1();
double analytical_solution_1(double time, double dot);

#endif //MFG_COMMON_H