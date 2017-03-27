#include <assert.h>

double *get_exact_solution(double *vec, int nx_1, double time, double hx, double a,
                           double (*analytical_solution)(double, double)) {
    assert(vec != NULL);
    assert(time > 0.);
    assert(nx_1 > 0);
    assert(hx > 0.);

    double *res = new double[nx_1];

    for (int i = 0; i < nx_1; i++)
        res[i] = analytical_solution(time, a + hx * i);

    return res;
}