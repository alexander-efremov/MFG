#ifndef MFG_COMMON_H
#define MFG_COMMON_H

#include "consts.h"

void get_exact_solution(double *res, int nx_1, double time, double hx, double a,double a_param,
                           double (*analytical_solution)(double, double, double));

double *solve_1();

#endif //MFG_COMMON_H