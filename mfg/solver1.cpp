#include <assert.h>
#include <stdio.h>
#include "common.h"
#include <string.h>
#include <malgo.h>
#include <utils.h>

/**
	 * n - число уравнений (строк матрицы)
	 * b - диагональ, лежащая под главной (нумеруется: [1;n-1])
	 * c - главная диагональ матрицы A (нумеруется: [0;n-1])
	 * d - диагональ, лежащая над главной (нумеруется: [0;n-2])
	 * f - правая часть (столбец)
	 * x - решение, массив x будет содержать ответ
	 */
inline void fill_b(double *b, int n) {
    b[0] = 0.;
    b[1] = 0.;
    for (int i = 2; i < n - 1; ++i)
        b[i] = 1. / (8. * TAU) - SIGMA_SQ / H_SQ;
    b[n - 1] = 0.;
}

inline void fill_c(double *c, int n) {
    double val = 7. / (8. * TAU) + SIGMA_SQ / H_SQ;

    c[0] = c[n - 1] = 0.;

    c[1] = c[n - 2] = val;

    for (int i = 2; i < n - 2; ++i)
        c[i] = 3. / (4. * TAU) + (2. * SIGMA_SQ) / H_SQ;
}

inline void fill_d(double *d, int n) {
    d[0] = 0.;
    for (int i = 1; i < n - 2; ++i)
        d[i] = 1. / (8. * TAU) - SIGMA_SQ / H_SQ;
    d[n - 1] = d[n - 2] = 0.;
}

inline double func_alpha(double a_coef, double t, double x) {
    return a_coef * t * x * (1. - x);
}

inline double get_f(double sigma_sq, double a_coef, double x, double t) {
    return (a_coef * x * x * (3. - 2. * x)) / 6.
           - (sigma_sq * a_coef * t * (1. - 2. * x)) / 2.
           + a_coef * a_coef * t * t * x * x * (1. - x) * (1. - x)
           + (a_coef * a_coef * t * t * x * x * (3. - 2. * x) * (1. - 2. * x)) / 6.;
}

double analytical_solution_1(double a, double t, double x) {
    return (a * t * x * x * (3. - 2. * x)) / 6.;
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
    alpha = func_alpha(A_COEF, time, x_left);
    x_left -= TAU * alpha;
    alpha = func_alpha(A_COEF, time, x_right);
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
    m_left = m_pr[I_left] * (xi_poh_left - x_left) / H + m_pr[I_left + 1] * (x_left - xi_moh_left) / H;

    // применим формулу трапеций - полусумма оснований умножить на высоту
    r += 0.5 * (m_left + m_pr[I_left + 1]) * (xi_poh_left - x_left);

    for (i = I_left + 1; i < I_right; ++i) {
        r += 0.5 * (m_pr[i] + m_pr[i + 1]) * H;
        printf("!!! INTERNAL POINTS  CALCULATION !!!\n");
    }

    // считаем интеграл по правой подчасти
    // вычислим значение функции в нашей правой точке по формуле 3.7
    m_right = m_pr[I_right] * (xi_poh_right - x_right) / H + m_pr[I_right + 1] * (x_right - xi_moh_right) / H;

    r += 0.5 * (m_right + m_pr[I_right]) * (x_right - xi_moh_right);

    return r;
}

void fill_rp(double *rp, double *m_pr, double time, int n) {
    rp[0] = 0.;
    rp[n - 1] = 0;

    for (int i = 1; i < n - 1; ++i)
        rp[i] = get_right_part_inner_points(i, m_pr, time);

    for (int i = 1; i < n - 1; ++i)
        rp[i] += get_f(SIGMA_SQ, A_COEF, (i - 0.5) * H, time);
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
    assert(A_COEF > 0.);
    assert(TIME_STEP_CNT >= 1);
    // (3.19)
//    printf("H * H = %e\n", H * H);
//    printf("8 * TAU * SIGMA_SQ = %e\n", 8. * TAU * SIGMA_SQ);
//    fflush(stdout);
    assert(H * H <= 8 * TAU * SIGMA_SQ);
}

/**
    * n - число уравнений (строк матрицы)
    * b - диагональ, лежащая под главной (нумеруется: [1;n-1])
    * c - главная диагональ матрицы A (нумеруется: [0;n-1])
    * d - диагональ, лежащая над главной (нумеруется: [0;n-2])
    * f - правая часть (столбец)
    * x - решение, массив x будет содержать ответ
    */
void print_thomas_arrays(double *b, double *c, double *d, int n) {
    printf("b\n");
    for (int i = 0; i < n; ++i) {
        printf("%e ", b[i]);
    }
    printf("\nc\n");
    for (int i = 0; i < n; ++i) {
        printf("%e ", c[i]);
    }
    printf("\nd\n");
    for (int i = 0; i < n; ++i) {
        printf("%e ", d[i]);
    }
    printf("\n");
}

void fill_arr_by_ex_sol(double *arr, int n, double time) {
    for (int i = 1; i < n - 1; ++i) arr[i] = analytical_solution_1(A_COEF, time, A + i * H - H_2);
    arr[0] = arr[1];
    arr[n - 1] = arr[n - 2];
}

double *solve_1(int n, double *exact_sol_to_fill) {
    assert_params();

    double *m = (double *) malloc(n * sizeof(double));
    double *m_pr = (double *) malloc(n * sizeof(double));
    double *rp = (double *) malloc(n * sizeof(double));
    double *ex_m = (double *) malloc(n * sizeof(double));
    double *err = (double *) malloc(n * sizeof(double));
    double *b = (double *) malloc(n * sizeof(double));
    double *c = (double *) malloc(n * sizeof(double));
    double *d = (double *) malloc(n * sizeof(double));
    fill_b(b, n);
    fill_c(c, n);
    fill_d(d, n);
    //print_thomas_arrays(b, c, d, n);

    for (int i = 0; i < n; ++i) m[i] = rp[i] = m_pr[i] = 0.;

    fill_arr_by_ex_sol(m_pr, n, 0.);

//    printf("M_PR\n");
//    print_matrix1(m_pr, 1, n);

    double *l1s = (double *) malloc(sizeof(double) * TIME_STEP_CNT);

    for (int tl = 1; tl <= TIME_STEP_CNT; ++tl) {
        fill_rp(rp, m_pr, TAU * tl, n);
//        printf("RP \n");
//        print_matrix1(rp, 1, n);
        thomas_algo_verzh_modified(n, b, c, d, rp, m);
        m[0] = m[1];
        m[n - 1] = m[n - 2];
        memcpy(m_pr, m, n * sizeof(double));
    }

//    printf("M DONE\n");
//    print_matrix1(m, 1, n);

    fill_arr_by_ex_sol(ex_m, n, TAU * TIME_STEP_CNT);
//    printf("EXACT SOL \n");
//    print_matrix1(ex_m, 1, n);

    fill_arr_diff(err, ex_m, m, n);
//    printf("ERR \n");
//    print_matrix1(err, 1, n);

    double l1 = get_l1_norm(H, n, err);
    printf("L1 NORM = %e\n", l1);

    fill_arr_by_ex_sol(exact_sol_to_fill, n, TAU * TIME_STEP_CNT);

    free(l1s);
    free(ex_m);
    free(m_pr);
    free(rp);
    free(err);
    free(b);
    free(c);
    free(d);

    return m;
}