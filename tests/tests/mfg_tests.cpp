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

void run_solver_1(unsigned int d) {
    A = 0.;
    A_COEF = 1.;
    B = 1.;
    NX = d;
    N_1 = NX + 1;
    H = (B - A) / NX;
    H_2 = 0.5 * H;
    H_SQ = H * H;
    SIGMA = 10.;
    SIGMA_SQ = SIGMA * SIGMA;
    TAU = 1.e-3;
    TIME_STEP_CNT = 1;

    print_params();

    double *density = solve_1();

    free(density);
}

TEST_F(MfgFixture, mfg_solver_1) {
    for (int i = 0; i < 1; ++i) {
        double d = 0;
        switch (i) {
            case 0:
                d = 10.;
                break;
            case 1:
                d = 100.;
                break;
            case 2:
                d = 200.;
                break;
            case 3:
                d = 400.;
                break;
            case 4:
                d = 800.;
                break;
            case 5:
                d = 1600.;
                break;
            default:
                return;
        }
        run_solver_1((unsigned int) d);
    }
}

