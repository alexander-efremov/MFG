#ifndef MFG_MALGO_H
#define MFG_MALGO_H

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <cfloat>

inline int is_diagonally_dominating_matrix(double *arr, int n, int m) {
    double diag_sum = 0.;
    double sum = 0.;

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            if (i == j) {
                diag_sum += fabs(arr[i * n + j]);
            } else {
                sum += fabs(arr[i * n + j]);
            }
        }
    }
    if (diag_sum >= sum) return 1;
    return 0;
}

inline void thomas_algo(int n, double *matrix, double *f, double *res) {
    assert(n > 2);
    assert(matrix != nullptr);
    assert(f != nullptr);
    assert(is_diagonally_dominating_matrix(matrix, n, n) == 1);

    double *alpha = (double *) malloc(n * sizeof(double));
    double *beta = (double *) malloc(n * sizeof(double));

    for (int i = 0; i < n; ++i) {
        res[i] = DBL_MIN;
        alpha[i] = DBL_MIN;
        beta[i] = DBL_MIN;
    }

    // step 1

    double y = matrix[0];
    alpha[0] = -matrix[1] / y;
    beta[0] = f[0] / y;

    int offset = n;
    for (int i = 1; i < n - 1; ++i) {

        double a_i = matrix[offset];
        double b_i = matrix[offset + 1];
        double c_i = matrix[offset + 2];

        y = b_i + a_i * alpha[i - 1];

        alpha[i] = -c_i / y;

        beta[i] = (f[i] - a_i * beta[i - 1]) / y;

        offset += n + 1;
    }

    double a_n = matrix[offset];
    double b_n = matrix[offset + 1];
    double alpha_n_minus_1 = alpha[n - 2];
    double beta_n_minus_1 = beta[n - 2];

    y = b_n + a_n * alpha_n_minus_1;
    beta[n - 1] = (f[n - 1] - a_n * beta_n_minus_1) / y;

    // step 2

    res[n - 1] = beta[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        res[i] = alpha[i] * res[i + 1] + beta[i];
    }

    free(alpha);
    free(beta);
}

/**
	 * n - число уравнений (строк матрицы)
	 * a - диагональ, лежащая под главной (нумеруется: [1;n-1])
	 * c - главная диагональ матрицы A (нумеруется: [0;n-1])
	 * b - диагональ, лежащая над главной (нумеруется: [0;n-2])
	 * f - правая часть (столбец)
	 * x - решение, массив x будет содержать ответ
	 */
inline void thomas_algo(int n, double *a, double *c, double *b, double *f, double *x) {
    double m;
    for (int i = 1; i < n; i++) {
        m = a[i] / c[i - 1];
        c[i] = c[i] - m * b[i - 1];
        f[i] = f[i] - m * f[i - 1];
    }

    x[n - 1] = f[n - 1] / c[n - 1];

    for (int i = n - 2; i >= 0; i--)
        x[i] = (f[i] - b[i] * x[i + 1]) / c[i];

}

inline void thomas_algo(double *x, int X, double *a, double *b, double *c) {
    /*
     solves Ax = v where A is a tridiagonal matrix consisting of vectors a, b, c
     x - initially contains the input vector v, and returns the solution x. indexed from 0 to X - 1 inclusive
     X - number of equations (length of vector x)
     a - subdiagonal (means it is the diagonal below the main diagonal), indexed from 1 to X - 1 inclusive
     b - the main diagonal, indexed from 0 to X - 1 inclusive
     c - superdiagonal (means it is the diagonal above the main diagonal), indexed from 0 to X - 2 inclusive

     Note: contents of input vector c will be modified, making this a one-time-use function (scratch space can be allocated instead for this purpose to make it reusable)
     Note 2: We don't check for diagonal dominance, etc.; this is not guaranteed stable
     */

    c[0] = c[0] / b[0];
    x[0] = x[0] / b[0];

    /* loop from 1 to X - 1 inclusive, performing the forward sweep */
    for (size_t ix = 1; ix < X; ix++) {
        const double m = 1.0f / (b[ix] - a[ix] * c[ix - 1]);
        c[ix] = c[ix] * m;
        x[ix] = (x[ix] - a[ix] * x[ix - 1]) * m;
    }

    /* loop from X - 2 to 0 inclusive (safely testing loop condition for an unsigned integer), to perform the back substitution */
    for (int ix = X - 2; ix >= 0; ix--)
        x[ix] -= c[ix] * x[ix + 1];
}

#endif //MFG_MALGO_H
