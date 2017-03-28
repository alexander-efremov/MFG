#include <assert.h>
#include <stdio.h>
#include "common.h"
#include <string.h>
#include <malgo.h>
#include <utils.h>

inline double get_lc_0() {
    return 7. / (8. * TAU) + SIGMA_SQ / (2. * H_SQ);
}

inline double get_lc_n() {
    return 7. / (8. * TAU) + SIGMA_SQ / (2. * H_SQ);
}

inline double get_lc_1() {
    return 1. / (8. * TAU) - SIGMA_SQ / (2. * H_SQ);
}

inline double get_lc_2() {
    return 3. / (4. * TAU) + SIGMA_SQ / H_SQ;
}

inline double get_lc_3() {
    return get_lc_1();
}

double *build_a(double *a, int n) {
    for (int i = 0; i < n * n; ++i) a[i] = 0.;

    a[0] = get_lc_0();
    a[1] = get_lc_3();

    int offset = n;
    for (int i = 1; i < n - 1; ++i) {
        a[offset] = get_lc_1();
        a[offset + 1] = get_lc_2();
        a[offset + 2] = get_lc_3();
        offset += n + 1;
    }

    a[n * n - 2] = get_lc_1();
    a[n * n - 1] = get_lc_n();

    return a;
}

inline double func_alpha(double a, double t, double x) {
    return a * t * x * (1 - x);
    //return ALPHA_COEF;
}

inline double get_rp_exact(double sigma_sq, double a, double x, double t) {
    return (a * x * x * (3. - 2. * x)) / 6.
           - (sigma_sq * a * t * (1. - 2. * x)) / 2.
           + a * a * t * t * x * x * (1. - x) * (1. - x)
           + (a * a * t * t * x * x * (3. - 2. * x) * (1. - 2. * x)) / 6.;
}

double analytical_solution_1(double a, double time, double x) {
    return a * time * (x * x / 6.) * (3 - 2 * x);
}

double get_right_part_inner_points(int I, double *m_pr, double time) {
    int j = I - 1;
    // moh - minus_one_half poh - plus_one_half
    double r = 0., m_left, m_right, xi_moh_left, xi_poh_left, xi_moh_right, xi_poh_right, x_left, x_right, alpha;
    int i, I_left, I_right;

    // определяем интервал в котором лежит наша точка
    x_left = A + j * H;
    x_right = A + (j + 1) * H;

    // опускаем траектории из этих точек
    alpha = func_alpha(ALPHA_COEF, time, x_left);
    x_left -= TAU * alpha;
    alpha = func_alpha(ALPHA_COEF, time, x_right);
    x_right -= TAU * alpha;
    if (x_left < A || x_left > B || x_right < A || x_right > B)
        printf("Time value %.8le! ERROR INDEX i=%d : x1=%.8le ** x2=%.8le\n ", time, j, x_left, x_right);

    I_left = (int) (((x_left - A) / H) + 0.5); // определяем индекс левого интервала линейности
    xi_moh_left = A + I_left * H - H_2; // левая граница левого интервала линейности
    xi_poh_left = A + I_left * H + H_2; // правая граница левого интервала линейности
    assert(x_left >= xi_moh_left && x_left <= xi_poh_left);

    I_right = (int) (((x_right - A) / H) + 0.5); // определяем индекс правого интервала линейности
    xi_moh_right = A + I_right * H - H_2; // левая граница правого интервала линейности
    xi_poh_right = A + I_right * H + H_2; // правая граница правого интервала линейности
    assert(x_right >= xi_moh_right && x_right <= xi_poh_right);

    // считаем интеграл по левой подчасти
    // вычислим значение функции в нашей левой точке по формуле 3.7
    if (x_left > H_2)
        m_left = m_pr[I_left] * (xi_poh_left - x_left) / H + m_pr[I_left + 1] * (x_left - xi_moh_left) / H;
    else
        m_left = m_pr[I_left];

    // применим формулу трапеций - полусумма оснований умножить на высоту
    r += 0.5 * (m_left + m_pr[I_left + 1]) * (A + (I_left + 1) * H - x_left);

    for (i = I_left + 1; i < I_right; ++i) {
        r += 0.5 * (m_pr[i] + m_pr[i + 1]) * H;
        printf("!!! INTERNAL POINTS  CALCULATION !!!");
    }

    // считаем интеграл по правой подчасти
    // вычислим значение функции в нашей правой точке по формуле 3.7
    if (x_right < B - H_2)
        m_right = m_pr[I_right] * (xi_poh_right - x_right) / H + m_pr[I_right + 1] * (x_right - xi_moh_right) / H;
    else
        m_right = m_pr[I_right + 1];

    r += 0.5 * (m_right + m_pr[I_right]) * (A + I_right * H - x_right);

    return r;
}

void fill_rp(double *rp, double *m_pr, double time) {
    // todo: ПЕРЕПИСАТЬ
    rp[0] = analytical_solution_1(ALPHA_COEF, time, A - H_2);
    rp[1] = analytical_solution_1(ALPHA_COEF, time, A);
    double val0 = get_rp_exact(SIGMA_SQ, ALPHA_COEF, (rp[0] + rp[1]) / 2., time);
    rp[0] += val0;
    rp[1] += val0;

    rp[N_1] = analytical_solution_1(ALPHA_COEF, time, A + N_1 * H);
    rp[N_1 + 1] = analytical_solution_1(ALPHA_COEF, time, A + (N_1 + 1) * H_2);
    double valN = get_rp_exact(SIGMA_SQ, ALPHA_COEF, (rp[N_1] + rp[N_1 + 1]) / 2., time);
    rp[N_1] += valN;
    rp[N_1 + 1] += valN;

    // inner points
    for (int i = 1; i < N_1; ++i) {
        double val = get_right_part_inner_points(i, m_pr, time);
        rp[i] = val;
//        if (rp[i] == 0.)
//            printf("rp error %d = %e\n", i, rp[i]);
    }

    // add F
    for (int i = 1; i < N_1; ++i) {
        double val = get_rp_exact(SIGMA_SQ, ALPHA_COEF, (rp[i - 1] + rp[i]) / 2., time);
        rp[i] = val;
    }
}

void fill_error(double *err, double *sol, int n, double time) {
    double *ex_sol = (double *) malloc(n * sizeof(double));
    for (int i = 0; i < n; ++i)
        ex_sol[i] = analytical_solution_1(ALPHA_COEF, time, A + i * H);

    print_matrix(ex_sol, 1, n);

    for (int i = 0; i < n; ++i)
        err[i] = ex_sol - sol;

    free(ex_sol);
}

void assert_params() {
    assert(H > 0.);
    assert(H_SQ == H * H);
    assert(H_2 == H / 2.);
    assert(SIGMA_SQ == SIGMA * SIGMA);
    assert(NX > 0);
    assert(N_1 == NX + 1);
    assert(A == 0.);
    assert(TAU > 0.);
    assert(ALPHA_COEF > 0.);
    assert(TIME_STEP_CNT >= 1);
    // (3.19)
//    printf("H * H = %e\n", H * H);
//    printf("8 * TAU * SIGMA_SQ = %e\n", 8 * TAU * SIGMA_SQ);
//    fflush(stdout);
    assert(H * H <= 8 * TAU * SIGMA_SQ);
}

double *solve_1() {
    assert_params();

    const unsigned int n = N_1 + 2;

    double *m = (double *) malloc(n * sizeof(double));
    double *m_pr = (double *) malloc(n * sizeof(double));
    double *rp = (double *) malloc(n * sizeof(double));
    double *a = (double *) malloc(N_1 * N_1 * sizeof(double));

    for (int i = 0; i < n; ++i) m[i] = rp[i] = m_pr[i] = 0.;

    m_pr[0] = analytical_solution_1(ALPHA_COEF, 0., A - H_2);
    for (int i = 1; i < N_1 + 1; ++i) m_pr[i] = analytical_solution_1(ALPHA_COEF, 0., A + i * H_2);
    m_pr[N_1 + 1] = analytical_solution_1(ALPHA_COEF, 0., B + H_2);
    print_matrix(m_pr, 1, n);

    for (int tl = 1; tl <= TIME_STEP_CNT; ++tl) {
        fill_rp(rp, m_pr, TAU * tl);
        //print_matrix(rp, 1, n);
        a = build_a(a, N_1);
        //print_matrix(a, N_1, N_1);
        thomas_algo(N_1, a, rp, m);
        memcpy(m_pr, m, n * sizeof(double));
    }
    print_matrix(m, 1, n);

    double *err = (double *) malloc(n * sizeof(double));
    fill_error(err, m, n, TIME_STEP_CNT * TAU);
    print_matrix(err, 1, n);
    free(err);

    free(m_pr);
    free(rp);
    free(a);

    return m;
}