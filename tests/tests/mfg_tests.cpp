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
    printf("\nNX = %dx\n", NX);
    printf("ALPHA = %le\n", ALPHA);
    printf("SIGMA = %le\n", SIGMA);
    printf("SIGMA_SQ = %le\n", SIGMA_SQ);
    printf("HX = %le\n", HX);
    printf("HX_SQ = %le\n", HX_SQ);
    printf("TAU = %le\n", TAU);
    printf("TIME_STEP_CNT = %d\n", TIME_STEP_CNT);
    fflush(stdout);
}

void run_solver_1(unsigned int d) {
    assert(d >= 50);

    A = 0.;
    B = 1.;
    SIGMA = 0.1;
    SIGMA_SQ = SIGMA * SIGMA;
    NX = d;
    NX_1 = NX + 1;
    HX = (B - A) / NX;
    HX_SQ = HX * HX;
    ALPHA = 1.;
    TAU = 1.e-3;
    TIME_STEP_CNT = 5;

    print_params();

    double *grid = (double*)malloc(NX_1* sizeof(double));
    double *gridPr = (double*)malloc(NX_1* sizeof(double));;

    double *density = solve_1(grid, gridPr);

    delete[] density;
    free(grid);
    free(gridPr);
}

TEST_F(MfgFixture, mfg_solver_1) {
    for (int i = 1; i < 2; ++i) {
        double d = 0;
        switch (i) {
            case 0:
                d = 50.;
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

