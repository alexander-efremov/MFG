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

double analytical_solution_1(double a, double time, double x) {
    return a * time * (x * x / 6.) * (3 - 2 * x);
}

double get_right_part_inner_points(int ind, double *m_pr, double time_value) {
    // moh - minus_one_half poh - plus_one_half
    double r = 0., m_left, m_right, xi_moh_left, xi_poh_left, xi_moh_right, xi_poh_right, x_left, x_right, u;
    int i, I_left, I_right;

    // определяем интервал в котором лежит наша точка
    x_left = A + (ind - 1) * H;
    x_right = A + ind * H;

    // опускаем траектории из этих точек
    u = func_alpha(ALPHA, time_value, x_left);
    x_left -= TAU * u;
    u = func_alpha(ALPHA, time_value, x_right);
    x_right -= TAU * u;
    if (x_left < A || x_left > B || x_right < A || x_right > B)
        printf("Time value %.8le! ERROR INDEX i=%d : x1=%.8le ** x2=%.8le\n ", time_value, ind, x_left, x_right);

    I_left = (int) (((x_left - A) / H) + 0.5); // определяем индекс левого интервала линейности
    xi_moh_left = A + I_left * H - 0.5 * H; // левая граница левого интервала линейности
    xi_poh_left = A + I_left * H + 0.5 * H; // правая граница левого интервала линейности

    I_right = (int) (((x_right - A) / H) + 0.5); // определяем индекс правого интервала линейности
    xi_moh_right = A + I_right * H - 0.5 * H; // левая граница правого интервала линейности
    xi_poh_right = A + I_right * H + 0.5 * H; // правая граница правого интервала линейности

    // считаем интеграл по левой подчасти
    // вычислим значение функции в нашей левой точке по формуле 3.7
    m_left = m_pr[I_left + 1] * (xi_poh_left - x_left) + m_pr[I_left + 2] * (x_left - xi_moh_left);
    // применим формулу трапеций - полусумма оснований умножить на высоту
    r += 0.5 * (m_left + m_pr[I_left + 1]) * (A + I_left * H - x_left);

    for (i = I_left + 1; i < I_right; ++i) r += 0.5 * (m_pr[i] + m_pr[i + 1]) * H;

    // считаем интеграл по правой подчасти
    // вычислим значение функции в нашей правой точке по формуле 3.7
    m_right = m_pr[I_right] * (xi_poh_right - x_right) + m_pr[I_right + 1] * (x_right - xi_moh_right);
    r += 0.5 * (m_right + m_pr[I_right]) * (A + I_right * H - x_right);

    return r;
}

void fill_rp(double *rp, double *m_pr, double time) {
    rp[0] = analytical_solution_1(ALPHA, time, A - 0.5 * H);
    rp[1] = analytical_solution_1(ALPHA, time, A);
    rp[N_1] = analytical_solution_1(ALPHA, time, A + N_1 * H);
    rp[N_1 + 1] = analytical_solution_1(ALPHA, time, A + (N_1 + 1) * H);

    // inner points
    for (int i = 2; i < N_1; ++i) {
        rp[i] = get_right_part_inner_points(i - 1, m_pr, time);
//        if (rp[i] == 0.)
//            printf("rp error %d = %e\n", i, rp[i]);
    }
}

double *solve_1() {
    double *m = (double *) malloc((N_1 + 2) * sizeof(double));
    double *m_pr = (double *) malloc((N_1 + 2) * sizeof(double));
    double *rp = (double *) malloc((N_1 + 2) * sizeof(double));
    double *a = (double *) malloc((N_1 + 2) * (N_1 + 2) * sizeof(double));

    for (int i = 0; i < N_1 + 2; ++i) m[i] = rp[i] = 0.;

    m_pr[0] = analytical_solution_1(ALPHA, 0., A - 0.5 * H);
    for (int i = 1; i < N_1 + 1; ++i) m_pr[i] = analytical_solution_1(ALPHA, 0., A + i * H);
    m_pr[N_1 + 1] = analytical_solution_1(ALPHA, 0., B + 0.5 * H);

    //print_matrix(m_pr, 1, N_1 + 2);
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