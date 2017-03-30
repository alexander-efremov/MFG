#ifndef MFG_MALGO_H
#define MFG_MALGO_H

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <cfloat>

/**
	 * n - число уравнений (строк матрицы)
	 * a - диагональ, лежащая под главной (нумеруется: [1;n-1])
	 * c - главная диагональ матрицы A (нумеруется: [0;n-1])
	 * b - диагональ, лежащая над главной (нумеруется: [0;n-2])
	 * f - правая часть (столбец)
	 * x - решение, массив x будет содержать ответ
	 */
inline void thomas_algo(int n, double *a, double *c, double *b, double *f, double *x) {
    assert(n > 2);
    assert(x != nullptr);
    assert(f != nullptr);
    assert(a != nullptr);
    assert(b != nullptr);
    assert(c != nullptr);

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

/**
	 * n - число уравнений (строк матрицы)
	 * b - диагональ, лежащая под главной (нумеруется: [1;n-1])
	 * c - главная диагональ матрицы A (нумеруется: [0;n-1])
	 * d - диагональ, лежащая над главной (нумеруется: [0;n-2])
	 * f - правая часть (столбец)
	 * x - решение, массив x будет содержать ответ
	 */
inline void thomas_algo_verzh(int n, double *b, double *c, double *d, double *r, double *x) {
    assert(n > 2);
    assert(x != nullptr);
    assert(r != nullptr);
    assert(b != nullptr);
    assert(c != nullptr);
    assert(d != nullptr);

    double *delta = (double *) malloc(n * sizeof(double));
    double *beta = (double *) malloc(n * sizeof(double));
    double *lambda = (double *) malloc(n * sizeof(double));

    for (int j = 0; j < n; ++j) delta[j] = beta[j] = lambda[j] = 0.;

    delta[0] = c[0];
    beta[0] = -d[0] / delta[0];
    lambda[0] = r[0] / delta[0];

    for (int i = 1; i < n; ++i) {
        delta[i] = c[i] + b[i] * beta[i - 1];
        beta[i] = -d[i] / delta[i];
        lambda[i] = (r[i] - b[i] * lambda[i - 1]) / delta[i];
    }

    x[n - 1] = lambda[n - 1];

    for (int i = n - 2; i >= 0; i--)
        x[i] = beta[i] * x[i + 1] + lambda[i];

    free(delta);
    free(beta);
    free(lambda);
}

/**
	 * n - число уравнений (строк матрицы)
	 * b - диагональ, лежащая под главной (нумеруется: [1;n-1])
	 * c - главная диагональ матрицы A (нумеруется: [0;n-1])
	 * d - диагональ, лежащая над главной (нумеруется: [0;n-2])
	 * f - правая часть (столбец)
	 * x - решение, массив x будет содержать ответ
	 */
inline void thomas_algo_verzh_modified(int n, double *b, double *c, double *d, double *r, double *x) {
    assert(n > 2);
    assert(x != nullptr);
    assert(r != nullptr);
    assert(b != nullptr);
    assert(c != nullptr);
    assert(d != nullptr);

    double *delta = (double *) malloc(n * sizeof(double));
    double *beta = (double *) malloc(n * sizeof(double));
    double *lambda = (double *) malloc(n * sizeof(double));

    for (int j = 0; j < n; ++j) delta[j] = beta[j] = lambda[j] = 0.;

    delta[1] = c[1];
    beta[1] = -d[1] / delta[1];
    lambda[1] = r[1] / delta[1];

    for (int i = 2; i < n-1; ++i) {
        delta[i] = c[i] + b[i] * beta[i - 1];
        beta[i] = -d[i] / delta[i];
        lambda[i] = (r[i] - b[i] * lambda[i - 1]) / delta[i];
    }

    x[n - 2] = lambda[n - 2];

    for (int i = n - 3; i >= 1; i--)
        x[i] = beta[i] * x[i + 1] + lambda[i];

    free(delta);
    free(beta);
    free(lambda);
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
