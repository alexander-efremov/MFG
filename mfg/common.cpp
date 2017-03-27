#include <assert.h>

void get_exact_solution(double *res, int nx_1, double time, double hx, double a, double a_param,
                           double (*analytical_solution)(double, double, double)) {
    assert(time > 0.);
    assert(nx_1 > 0);
    assert(hx > 0.);

    for (int i = 0; i < nx_1; i++)
        res[i] = (*analytical_solution)(a_param, time, a + hx * i);
}