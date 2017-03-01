#include <utils.h>
#include "gtest/gtest.h"

class Fixture : public ::testing::Test {
protected:
    virtual void TearDown() {

    }

    virtual void SetUp() {
    }

public:
    Fixture() : Test() {
    }
};

void print_params() {

    fflush(stdout);
}

TEST_F(Fixture, thomas_algo_test) {

}

