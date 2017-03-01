#ifndef MFG_MALGO_H
#define MFG_MALGO_H

#include <stdlib.h>
#include <assert.h>
#include <cfloat>

inline double *thomas_algo(int n, double *arr, double *f) {
    assert(n > 2);
    assert(arr != nullptr);
    assert(f != nullptr);

    double *res = (double *) malloc(n * sizeof(double));
    double *alpha = (double *) malloc(n * sizeof(double));
    double *beta = (double *) malloc(n * sizeof(double));

    for (int i = 0; i < n; ++i) {
        res[i] = DBL_MIN;
        alpha[i] = DBL_MIN;
        beta[i] = DBL_MIN;
    }

    // step 1

    double y = arr[0];
    alpha[0] = -arr[1] / y;
    beta[0] = f[0] / y;

    int offset = n;
    for (int i = 1; i < n - 1; ++i) {

        double a_i = arr[offset];
        double b_i = arr[offset + 1];
        double c_i = arr[offset + 2];

        y = b_i + a_i * alpha[i - 1];

        alpha[i] = -c_i / y;

        beta[i] = (f[i] - a_i * beta[i - 1]) / y;

        offset += n + 1;
    }

    double b_n = arr[n - 1];
    double a_n = arr[n - 2];

    double alpha_n_minus_1 = alpha[n - 2];
    double beta_n_minus_1 = beta[n - 2];

    y = b_n + a_n * alpha_n_minus_1;
    beta[n - 1] = (f[n - 1] - a_n * beta_n_minus_1) / y;

    // step 2

    res[n - 1] = beta[n - 1];
    for (int i = n - 2; i < 0; ++i) {
        res[i] = alpha[i] * res[i + 1] + beta[i];
    }

    free(alpha);
    free(beta);

    return res;
}

#endif //MFG_MALGO_H
