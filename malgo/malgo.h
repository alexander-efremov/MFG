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
    //assert(is_diagonally_dominating_matrix(matrix, n, n) == 1);

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

#endif //MFG_MALGO_H
