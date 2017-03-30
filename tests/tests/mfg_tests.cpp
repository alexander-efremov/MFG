#include <consts.h>
#include <common.h>
#include <utils.h>
#include "gtest/gtest.h"

class MfgFixture : public ::testing::Test {
protected:
    virtual void TearDown() {
    }

    virtual void SetUp() {
    }

public:
    MfgFixture() : Test() {
    }
};

void print_params() {
    printf("\nNX = %d\n", NX);
    printf("A_COEF = %le\n", A_COEF);
    printf("SIGMA = %le\n", SIGMA);
    printf("SIGMA_SQ = %le\n", SIGMA_SQ);
    printf("H = %le\n", H);
    printf("H_SQ = %le\n", H_SQ);
    printf("TAU = %le\n", TAU);
    printf("TIME_STEP_CNT = %d\n", TIME_STEP_CNT);
//    printf("%f ", A - H_2);
//    for (int i = 1; i < N_1; ++i) {
//        printf("%f ", i * H);
//    }
//    printf("%f\n", B + H_2);
    fflush(stdout);
}

double run_solver_1(int nx, double tau, int time_step_cnt, double sigma) {
    A = 0.;
    A_COEF = 1.;
    B = 1.;
    SIGMA = sigma;
    SIGMA_SQ = SIGMA * SIGMA;


    NX = nx;
    N_1 = NX + 1;
    H = (B - A) / NX;
    H_2 = 0.5 * H;
    H_SQ = H * H;

    TAU = tau;
    TIME_STEP_CNT = time_step_cnt;

    if (sigma == 1.)
        assert(tau >= H_SQ / 8.);

    print_params();

    int n = N_1 + 1;
    double *exact_m = (double *) malloc(n * sizeof(double));
    double *err = (double *) malloc(n * sizeof(double));
    double *m = solve_1(n, exact_m);

    for (int i = 0; i < n; ++i)
        err[i] = exact_m[i] - m[i];

    double l1 = get_l1_norm(H, n, err);

    free(m);
    free(exact_m);
    free(err);

    return l1;
}

TEST_F(MfgFixture, mfg_solver_1) {
    int exp_cnt = 6;
    double *l1_arr = (double *) malloc(exp_cnt * sizeof(double));
    for (int i = 0; i < exp_cnt; ++i) {
        printf("\n\n =============== EXPERIMENT %d ================ \n\n", i + 1);
        int nx = 0;
        double tau = 1e-3;
        double sigma = 1.;
        int tsc = 10;
        double mult = 1.;
        switch (i) {
            case 0:
                nx = 50;
                tau = tau / mult;
                tsc = 10 * (int) mult;
                break;
            case 1:
                mult = 4.;
                nx = 100;
                tau = tau / mult;
                tsc = 10 * (int) mult;
                break;
            case 2:
                mult = 16.;
                nx = 200;
                tau = tau / mult;
                tsc = 10 * (int) mult;
                break;
            case 3:
                mult = 32.;
                nx = 400;
                tau = tau / mult;
                tsc = 10 * (int) mult;
                break;
            case 4:
                mult = 64.;
                nx = 800;
                tau = tau / mult;
                tsc = 10 * (int) mult;
                break;
            case 5:
                mult = 128.;
                nx = 1600;
                tau = tau / mult;
                tsc = 10 * (int) mult;
                break;
            default:
                return;
        }
        printf("TAU * TIME_STEP_COUNT = %e", tau * tsc);
        assert(tau * tsc == 0.01);//0.01 sec
        double l1 = run_solver_1(nx, tau, tsc, sigma);
        l1_arr[i] = l1;

    }

    double l1_max = -10000000.;
    for (int j = 0; j < exp_cnt; ++j) {
        if (l1_arr[j] > l1_max)
            l1_max = l1_arr[j];
    }
    printf("\n\nL1 MAX ON EXPERIMENTS: %e\n", l1_max);

    free(l1_arr);
}

