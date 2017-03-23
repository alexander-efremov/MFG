#include <assert.h>
#include <stdio.h>
#include "common.h"
#include <string.h>
#include <malgo.h>
#include <utils.h>

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

inline double func_alpha(double a, double t, double x) {
    return a * t * x * (1 - x);
    //return ALPHA;
}

double get_right_part_inner_points(int ii, double *m_pr, double time_value) {

    double x_left = A + ii * HX;
    double x_right = A + (ii + 1) * HX;

    double u = func_alpha(ALPHA, time_value, x_left);
    x_left -= TAU * u;
    u = func_alpha(ALPHA, time_value, x_right);
    x_right -= TAU * u;
    if (x_left <= A || x_left >= B || x_right <= A || x_right >= B)
        printf("Time value %.8le! ERROR INDEX i=%d : x1=%.8le ** x2=%.8le\n ", time_value, ii, x_left, x_right);

    double r = -1.;
    double m_left;
    double m_right;

    int int_left = (int) (((x_left - A) / HX) + 0.5);
    int int_right = (int) (((x_right - A) / HX) + 0.5);

    double xi_plus_one_half = A + int_left * HX + HX / 2.;
    double xi_minus_one_half = A + int_left * HX - HX / 2.;
    m_left = m_pr[int_left] * (xi_plus_one_half - x_left) + m_pr[int_left] * (xi_minus_one_half - x_left);
    double val_left = 0.5 * (m_left + m_pr[int_left + 1]) * (A + int_left * HX - x_left);
    r += val_left;

    for (int i = int_left + 1; i < int_right; ++i) {
        double val_middle = 0.5 * (m_pr[i] + m_pr[i + 1]) * HX;
        r += val_middle;
    }

    xi_plus_one_half = A + int_right * HX + HX / 2.;
    xi_minus_one_half = A + int_right * HX - HX / 2.;
    m_right = m_pr[int_right] * (xi_plus_one_half - x_right) + m_pr[int_right] * (xi_minus_one_half - x_right);
    double val_right = 0.5 * (m_right + m_pr[int_right + 1]) * (A + int_right * HX - x_right);
    r += val_right;

    return r;
}

void fill_rp(double *rp, double *m_pr, double time) {
    rp[0] = analytical_solution_1(ALPHA, 0., A + 0 * HX);
    rp[1] = analytical_solution_1(ALPHA, 0., A + 0.5 * HX);
    rp[N_1] = analytical_solution_1(ALPHA, 0., A + N_1 * HX);;
    rp[N_1 + 1] = analytical_solution_1(ALPHA, 0., A + (N_1 + 1) * HX);;

    // inner points
    for (int i = 1; i < N_1 + 1; ++i) {
        rp[i] = get_right_part_inner_points(i, m_pr, time);
        if (rp[i] == -1.)
            printf("rp error %d = %e\n", i, rp[i]);
    }
}

double analytical_solution_1(double a, double time, double x) {
    return a * time * (x * x / 6.) * (3 - 2 * x);
}

double *solve_1() {
    double *m = (double *) malloc((N_1 + 2) * sizeof(double));
    double *m_pr = (double *) malloc((N_1 + 2) * sizeof(double));
    double *rp = (double *) malloc(N_1 * sizeof(double));
    double *a = (double *) malloc(N_1 * N_1 * sizeof(double));

    for (int i = 0; i < NX + 2; ++i) {
        m[i] = 0.;
        m_pr[i] = 0.;
    }

    for (int i = 0; i < N_1 + 2; ++i) {
        m_pr[i] = analytical_solution_1(ALPHA, 0., A + i * HX);
    }

    for (int tl = 1; tl <= TIME_STEP_CNT; ++tl) {
        fill_rp(rp, m_pr, TAU * tl);

        print_matrix(rp, 1, N_1 + 2);

        a = build_a(N_1);

        thomas_algo(N_1, a, rp, m);

        memcpy(m_pr, m, (N_1 + 2) * sizeof(double));
    }

    free(m_pr);
    free(rp);
    free(a);

    return m;
}




