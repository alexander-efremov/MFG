#include <malgo.h>
#include <utils.h>
#include "gtest/gtest.h"

class MiscFixture : public ::testing::Test {
protected:
    virtual void TearDown() {
    }

    virtual void SetUp() {
    }

public:
    MiscFixture() : Test() {
    }
};

//inline double get_lc_0() {
//    return 3;
//}
//
//inline double get_lc_n() {
//    return 1;
//}
//
//inline double get_lc_1() {
//    return 1;
//}
//
//inline double get_lc_2() {
//    return 2;
//}
//
//inline double get_lc_3() {
//    return 1;
//}
//
//double *build_a1(int n) {
//    double *r = (double *) malloc(n * n * sizeof(double));
//    for (int i = 0; i < n * n; ++i) {
//        r[i] = 0.;
//    }
//
//    r[0] = get_lc_0();
//    r[1] = get_lc_2();
//
//    int offset = n;
//    for (int i = 1; i < n - 1; ++i) {
//        r[offset] = get_lc_1();
//        r[offset + 1] = get_lc_2();
//        r[offset + 2] = get_lc_3();
//        offset += n + 1;
//    }
//
//    r[n * n - 2] = get_lc_n();
//    r[n * n - 1] = get_lc_2();
//    return r;
//}

TEST_F(MiscFixture, build_a1) {
//    const int n = 50;
//    double *r = build_a1(n);
//    printf("\n");
//    printf("\n");
//    print_matrix(r, n, n, 3);
//    ASSERT_TRUE(is_diagonally_dominating_matrix(r, n, n) == 1);
//
//    printf("\n");
//    free(r);
}


