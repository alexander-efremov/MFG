#include <utils.h>
#include <common.h>
#include "gtest/gtest.h"

class FemFixture : public ::testing::Test {
protected:
    virtual void TearDown() {

    }

    virtual void SetUp() {
    }

public:
    FemFixture() : Test() {
    }
};

void print_params() {

    fflush(stdout);
}

void print_params_1() {

    fflush(stdout);
}

void init_boundary_arrays_and_cp(int nx, int ny) {

}

void run_solver_1(unsigned int d) {
    assert(d >= 50);

}

// убрано притягивание сетки
// убрана рекурсия
// деревья вместо массивов?
//
TEST_F(FemFixture, test1) {
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

