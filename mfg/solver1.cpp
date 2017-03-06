#include <assert.h>
#include <math.h>
#include "common.h"

inline double get_lc_1() {
    return 1. / (8. * TAU) - SIGMA_SQ / (2. * HX_SQ);
}

inline double get_lc_2() {
    return 3. / (4. * TAU) + SIGMA_SQ / HX_SQ;
}

inline double get_lc_3() {
    return get_lc_1();
}

inline double func_alpha(double t, double x) {
    return ALPHA;
}

double *solve_1(int *grid, int *gridPr) {
    assert(grid != NULL);
    assert(gridPr != NULL);
    assert(HX_SQ = sqrt(HX));
    assert(SIGMA_SQ = sqrt(SIGMA));

    return nullptr;
}

double *calc_exact_1(int *grid, double t, int nx_1, double hx) {
    assert(grid != NULL);
    assert(t > 0.);
    assert(nx_1 > 0);
    assert(hx > 0.);
    return nullptr;
}