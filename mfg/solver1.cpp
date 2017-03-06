#include <assert.h>
#include <math.h>
#include "common.h"
#include "malgo.h"

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

double analytical_solution_1(double time, double dot) {
    return sin(time * dot);
}

double *solve_1(double *vec, double *vecPr) {
    assert(vec != NULL);
    assert(vecPr != NULL);
    assert(HX_SQ = sqrt(HX));
    assert(SIGMA_SQ = sqrt(SIGMA));

    for (int j = 0; j < NX_1; ++j) {
        vecPr[j] = analytical_solution_1(0, vecPr[j]);
    }
    // corner point i = 0;

    // corner point i = N;

    // inner points
    for (int i = 1; i < NX; ++i) {

    }
    return nullptr;
}



