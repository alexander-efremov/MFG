#include <malgo.h>
#include <utils.h>
#include "gtest/gtest.h"

class MalgoFixture : public ::testing::Test {
public:
    MalgoFixture() : Test() {
    }
};

TEST_F(MalgoFixture, thomas_algo_test3) {
    const int n = 3;

    double *a = (double *) malloc(n * sizeof(double));
    double *b = (double *) malloc(n * sizeof(double));
    double *c = (double *) malloc(n * sizeof(double));
    double *rp = (double *) malloc(n * sizeof(double));
    double *m = (double *) malloc(n * sizeof(double));
    for (int i = 0; i < n; ++i)
        a[i] = b[i] = c[i] = rp[i] = m[i] = 0.;

    //под
    a[0] = 0;
    a[1] = 5;
    a[2] = 1;

    //главная
    b[0] = 2;
    b[1] = 4;
    b[2] = -3;

    //над
    c[0] = -1;
    c[1] = 2;
    c[2] = 0;

    rp[0] = 3;
    rp[1] = 6;
    rp[2] = 2;

    printf("\n");
    print_matrix(rp, 1, n);

    thomas_algo(rp, n, a, b, c);

    print_matrix(rp, 1, n);

    ASSERT_NEAR(rp[0], 1.49, 1e-2);
    ASSERT_NEAR(rp[1], -0.02, 1e-2);
    ASSERT_NEAR(rp[2], -0.68, 1e-2);

    free(rp);
    free(a);
    free(b);
    free(c);
}

TEST_F(MalgoFixture, thomas_algo_test4) {
    const int n = 3;

    double *a = (double *) malloc(n * sizeof(double));
    double *b = (double *) malloc(n * sizeof(double));
    double *c = (double *) malloc(n * sizeof(double));
    double *rp = (double *) malloc(n * sizeof(double));
    double *m = (double *) malloc(n * sizeof(double));
    for (int i = 0; i < n; ++i) {
        a[i] = b[i] = c[i] = rp[i] = m[i] = 0.;
    }

    //под
    a[0] = 0;
    a[1] = 1;
    a[2] = 1;

    //главная
    b[0] = 1;
    b[1] = 2;
    b[2] = 3;

    //над
    c[0] = 1;
    c[1] = 1;
    c[2] = 0;

    rp[0] = 2;
    rp[1] = 4;
    rp[2] = 4;

    printf("\n");
    print_matrix(rp, 1, n);

    thomas_algo(rp, n, a, b, c);

    print_matrix(rp, 1, n);

    ASSERT_NEAR(rp[0], 1., 1e-2);
    ASSERT_NEAR(rp[1], 1., 1e-2);
    ASSERT_NEAR(rp[2], 1., 1e-2);

    free(rp);
    free(a);
    free(b);
    free(c);
}

TEST_F(MalgoFixture, thomas_algo_test5) {
    const int n = 3;

    double *a = (double *) malloc(n * sizeof(double));
    double *b = (double *) malloc(n * sizeof(double));
    double *c = (double *) malloc(n * sizeof(double));
    double *rp = (double *) malloc(n * sizeof(double));
    double *m = (double *) malloc(n * sizeof(double));
    for (int i = 0; i < n; ++i) {
        a[i] = b[i] = c[i] = rp[i] = m[i] = 0.;
    }

    //под
    a[0] = 0;
    a[1] = 5;
    a[2] = 1;

    //главная
    c[0] = 2;
    c[1] = 4;
    c[2] = -3;

    //над
    b[0] = -1;
    b[1] = 2;
    b[2] = 0;

    rp[0] = 3;
    rp[1] = 6;
    rp[2] = 2;

    printf("\n");
    print_matrix(rp, 1, n);
/**
	 * n - число уравнений (строк матрицы)
	 * a - диагональ, лежащая под главной (нумеруется: [1;n-1])
	 * c - главная диагональ матрицы A (нумеруется: [0;n-1])
	 * b - диагональ, лежащая над главной (нумеруется: [0;n-2])
	 * f - правая часть (столбец)
	 * x - решение, массив x будет содержать ответ
	 */
    thomas_algo(n, a, c, b, rp, m);

    print_matrix(m, 1, n);

    ASSERT_NEAR(m[0], 1.49, 1e-2);
    ASSERT_NEAR(m[1], -0.02, 1e-2);
    ASSERT_NEAR(m[2], -0.68, 1e-2);

    free(rp);
    free(a);
    free(b);
    free(c);
}

TEST_F(MalgoFixture, thomas_algo_test6) {
    const int n = 3;

    double *a = (double *) malloc(n * sizeof(double));
    double *b = (double *) malloc(n * sizeof(double));
    double *c = (double *) malloc(n * sizeof(double));
    double *rp = (double *) malloc(n * sizeof(double));
    double *m = (double *) malloc(n * sizeof(double));
    for (int i = 0; i < n; ++i)
        a[i] = b[i] = c[i] = rp[i] = m[i] = 0.;

    //под
    a[0] = 0;
    a[1] = 1;
    a[2] = 1;

    //главная
    c[0] = 1;
    c[1] = 2;
    c[2] = 3;

    //над
    b[0] = 1;
    b[1] = 1;
    b[2] = 0;

    rp[0] = 2;
    rp[1] = 4;
    rp[2] = 4;

    printf("\n");
    print_matrix(rp, 1, n);

    thomas_algo(n, a, c, b, rp, m);

    print_matrix(m, 1, n);

    ASSERT_NEAR(m[0], 1., 1e-2);
    ASSERT_NEAR(m[1], 1., 1e-2);
    ASSERT_NEAR(m[2], 1., 1e-2);

    free(rp);
    free(a);
    free(b);
    free(c);
}

TEST_F(MalgoFixture, thomas_algo1_test1) {
    const int n = 5;

    double *b = (double *) malloc(n * sizeof(double));
    double *c = (double *) malloc(n * sizeof(double));
    double *d = (double *) malloc(n * sizeof(double));
    double *r = (double *) malloc(n * sizeof(double));
    double *x = (double *) malloc(n * sizeof(double));
    for (int i = 0; i < n; ++i)
        b[i] = b[i] = c[i] = d[i] = x[i] = 0.;

    //под
    b[0] = 0.;
    b[1] = 2.;
    b[2] = 4.;
    b[3] = 4.;
    b[4] = 2.;

    //главная
    c[0] = 2.;
    c[1] = 9.;
    c[2] = 17.;
    c[3] = 15.;
    c[4] = 3.;

    //над
    d[0] = 1.;
    d[1] = 2.;
    d[2] = -4.;
    d[3] = -8.;
    d[4] = 0.;

    r[0] = -10.;
    r[1] = -26.;
    r[2] = -16.;
    r[3] = -2.;
    r[4] = 16.;

    printf("\n");
    print_matrix1(b, 1, n);
    printf("\n");
    print_matrix1(c, 1, n);
    printf("\n");
    print_matrix1(d, 1, n);
    printf("\n");
    print_matrix1(r, 1, n);

    thomas_algo_verzh(n, b, c, d, r, x);

    printf("\n");
    print_matrix(x, 1, n);

    ASSERT_NEAR(x[0], -4., 1e-2);
    ASSERT_NEAR(x[1], -2., 1e-2);
    ASSERT_NEAR(x[2], 0., 1e-2);
    ASSERT_NEAR(x[3], 2., 1e-2);
    ASSERT_NEAR(x[4], 4., 1e-2);

    free(r);
    free(b);
    free(c);
    free(d);
    free(x);
}

TEST_F(MalgoFixture, thomas_algo1_test2) {
    const int n = 3;

    double *b = (double *) malloc(n * sizeof(double));
    double *c = (double *) malloc(n * sizeof(double));
    double *d = (double *) malloc(n * sizeof(double));
    double *r = (double *) malloc(n * sizeof(double));
    double *x = (double *) malloc(n * sizeof(double));
    for (int i = 0; i < n; ++i)
        b[i] = b[i] = c[i] = d[i] = x[i] = 0.;

    //под
    b[0] = 0.;
    b[1] = 1.;
    b[2] = 1.;

    //главная
    c[0] = 1.;
    c[1] = 2.;
    c[2] = 3.;

    //над
    d[0] = 1.;
    d[1] = 1.;
    d[2] = 0.;

    r[0] = 2.;
    r[1] = 4.;
    r[2] = 4.;

    printf("\n");
    print_matrix1(b, 1, n);
    printf("\n");
    print_matrix1(c, 1, n);
    printf("\n");
    print_matrix1(d, 1, n);
    printf("\n");
    print_matrix1(r, 1, n);

    thomas_algo_verzh(n, b, c, d, r, x);

    printf("\n");
    print_matrix(x, 1, n);

    ASSERT_NEAR(x[0], 1., 1e-2);
    ASSERT_NEAR(x[1], 1., 1e-2);
    ASSERT_NEAR(x[2], 1., 1e-2);

    free(r);
    free(b);
    free(c);
    free(d);
    free(x);
}