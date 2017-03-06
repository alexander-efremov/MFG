#include <assert.h>
#include <stdio.h>
#include "common.h"
#include <math.h>
#include <string.h>
#include <malgo.h>

inline double get_lc_0() {
    return (7. / (8. * TAU)) - (SIGMA_SQ / (2. * HX_SQ));
}

inline double get_lc_n() {
    return (7. / (8. * TAU)) - (SIGMA_SQ / (2. * HX_SQ));
}

inline double get_lc_1() {
    return 1. / (8. * TAU) - SIGMA_SQ / (2. * HX_SQ);
}

inline double get_lc_2() {
    return 3. / (4. * TAU) + SIGMA_SQ / HX_SQ;
}

inline double get_lc_3() {
    return get_lc_1();
}

double *build_a(int n) {
    double *r = (double *) malloc(n * n * sizeof(double));
    for (int i = 0; i < n * n; ++i) {
        r[i] = 0.;
    }

    r[0] = get_lc_0();
    r[1] = get_lc_2();

    int offset = n;
    for (int i = 1; i < n - 1; ++i) {
        r[offset] = get_lc_1();
        r[offset + 1] = get_lc_2();
        r[offset + 2] = get_lc_3();
        offset += n + 1;
    }

    r[n * n - 2] = get_lc_n();
    r[n * n - 1] = get_lc_2();
    return r;
}

inline double func_alpha(double t, double x) {
    return ALPHA;
}

double get_right_part_inner_points(int ii, double *m_pr, double time_value) {

    double x_hat_left = A + ii * HX;
    double x_hat_right = A + (ii + 1) * HX;

    double u = func_alpha(time_value, x_hat_left);
    x_hat_left = x_hat_left - TAU * u;
    u = func_alpha(time_value, x_hat_right);
    x_hat_right = x_hat_right - TAU * u;
    if (x_hat_left <= A || x_hat_left >= B || x_hat_right <= A || x_hat_right >= B)
        printf("Time value %.8le! ERROR INDEX i=%d : x1=%.8le ** x2=%.8le\n ", time_value, ii, x_hat_left, x_hat_right);

    int sq_i_left = (int) ((x_hat_left - A) / HX);
    int sq_i_right = (int) ((x_hat_right - A) / HX);

    double r = -1.;
    if (sq_i_right - 1 == sq_i_left) // intervals are joint
    {
        double xi_plus_one_half = A + sq_i_left * HX + HX / 2.;
        double xi_minus_one_half = A + sq_i_left * HX - HX / 2.;

        int idx1 = (int) ((xi_plus_one_half - x_hat_left) / HX);
        int idx2 = (int) ((x_hat_left - xi_minus_one_half) / HX);
        r = idx1 + idx2;

        xi_plus_one_half = A + sq_i_right * HX + HX / 2.;
        xi_minus_one_half = A + sq_i_right * HX - HX / 2.;

        idx1 = (int) ((xi_plus_one_half - x_hat_right) / HX);
        idx2 = (int) ((x_hat_right - xi_minus_one_half) / HX);
        r += idx1 + idx2;

    } else {
        int diff = sq_i_right - sq_i_left;

        while (diff > 0) {
            diff--;
            // TODO
            printf("TODO: INTERVALS ISN'T JOINDED");
        }

    }


    return r;
}

double *get_rp(double *rp, double *m_pr, double time) {
    // i = 0;
    // TODO
    // i = N;
    // TODO
    // inner points
    for (int i = 2; i < NX * 2 + 1; ++i)
        rp[i] = get_right_part_inner_points(i, m_pr, time);
}


double analytical_solution_1(double time, double dot) {
    return 1. + sin(time * dot);
}

double *solve_1() {
    double *m = (double *) malloc((NX_1 * 2 + 2) * sizeof(double));
    double *m_pr = (double *) malloc((NX_1 * 2 + 2) * sizeof(double));
    double *rp = (double *) malloc(NX_1 * sizeof(double));
    double *a = (double *) malloc(NX_1 * NX_1 * sizeof(double));

    for (int i = 0; i < NX * 2 + 2; ++i) {
        m[i] = 0.;
        m_pr[i] = 0.;
    }

    for (int i = 1; i < NX * 2 + 1; ++i) {
        m_pr[i] = analytical_solution_1(0., A + i * HX);
    }

    for (int tl = 1; tl <= TIME_STEP_CNT; tl++) {
        rp = get_rp(rp, m_pr, TAU * tl);
        a = build_a(NX_1);

        thomas_algo(NX_1, a, rp, m);

        memcpy(m_pr, m, NX_1 * sizeof(double));
    }

    free(m_pr);
    free(rp);
    free(a);

    return m;
}




