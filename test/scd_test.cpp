#include <gtest/gtest.h>
#include <cmath>
#include <fstream>
#include "../include/scd/scd.h"
#include "data_utilities.h"
#include <iostream>

GTEST_TEST(Test2D, CreateShape){
    Eigen::Rotation2D<double> Xr(-M_PI/5);
    Eigen::Translation<double,2> Xt(1,-1);
    Eigen::Transform<double,2,Eigen::Isometry> X = Xt*Xr;
    scd::Superquadric2 SQ1(1,2,0.25);
    SQ1.X = X;

    Xr.angle() = -M_PI/3;
    Xt.x() = -3;
    Xt.y() = 3;
    X = Xt*Xr;
    scd::Superquadric2 E2(0.5,2.0,1.0);
    E2.X = X;

    scd::Request2 request;
    scd::Result2 result;
    scd::collide2(SQ1, E2, request, result);

    Eigen::Transform<double,2,Eigen::Isometry> X_expect;
    X_expect.matrix() <<  0.5000,    0.8660,   -0.3587,
                         -0.8660,    0.5000,    0.3587,
                          0,         0,         1.0000;

    ASSERT_LT(( X_expect.matrix() - result.Ec.X.matrix() ).norm(), 1E-4);
    EXPECT_NEAR(result.theta,2.726233503871522, 1E-6);
    ASSERT_FALSE(result.collide);
}
GTEST_TEST(TestCollision, Data){
    std::string file_path = "../test/data/cube_ellipsoid.txt";
    scd::Superquadric3 SQ1;
    scd::Ellipsoid3 E2;
    scd::Ellipsoid3 E2c_expect;
    scd::Result result_expect;
    uint n_line = 0;
    scd::test::ReadTestFile(file_path, n_line, SQ1, E2, E2c_expect, result_expect);

    scd::Request request;
    scd::Result result;
    scd::collide(SQ1, E2, request, result);

    ASSERT_EQ(result.collide, result_expect.collide);
    EXPECT_NEAR(result.eta, result_expect.eta, 1E-6);
    EXPECT_NEAR(result.omega, result_expect.omega, 1E-6);
}
GTEST_TEST(TestCollision, EllipsoidEllipsoidBatch){
    std::string file_path = "../test/data/ellipsoid_ellipsoid_batch.txt";
    uint n_line = 0;
    std::cout << "Begin batch test...\n";
    for (int n_test = 0; n_test < 1000; n_test++) {
        scd::Superquadric3 SQ1;
        scd::Ellipsoid3 E2;
        scd::Ellipsoid3 E2c_expect;
        scd::Result result_expect;
        scd::test::ReadTestFile(file_path, n_line, SQ1, E2, E2c_expect, result_expect);
        scd::Request request;
        scd::Result result;
        scd::collide(SQ1, E2, request, result);

        ASSERT_EQ(result.collide, result_expect.collide);
        // TODO (jay) Add angle wrapping to allow comparison
        n_line += 25;
    }
}

GTEST_TEST(TestCollision, CubeEllipsoid){
    scd::Superquadric3 SQ1(1.0,2.0,3.0,0.25,0.25);
    auto X_SQ1_R = Eigen::AngleAxis<double>(M_PI/3, Eigen::Vector3d::UnitX())
                 * Eigen::AngleAxis<double>(2*M_PI/3, Eigen::Vector3d::UnitY());
    auto X_SQ1_t = Eigen::Translation<double,3>(-1.0,1.0,2.0);
    SQ1.X = X_SQ1_t*X_SQ1_R;

    scd::Ellipsoid3 E2(1.0,2.0,0.5);
    auto X_E2_R = Eigen::AngleAxis<double>(7*M_PI/6, Eigen::Vector3d::UnitX())
                * Eigen::AngleAxis<double>(-4*M_PI/6, Eigen::Vector3d::UnitY())
                * Eigen::AngleAxis<double>(3*M_PI/5, Eigen::Vector3d::UnitZ());
    auto X_E2_t = Eigen::Translation<double,3>(9.0,3.0,6.0);
    E2.X = X_E2_t*X_E2_R;

    scd::Request request;
    scd::Result result;
    scd::collide(SQ1, E2, request, result);

    Eigen::Transform<double,3,Eigen::Isometry> X;
    X.matrix() <<      0.1545,    0.4755,   -0.8660,    2.9909,
                      -0.9574,   -0.1442,   -0.2500,    1.7982,
                      -0.2438,    0.8678,    0.4330,    3.5964,
                       0,         0,         0,         1.0000;

    // Data was copied from MATLAB printout, keeping only 4 decimals. Tolerances reflect this.
    ASSERT_LT(( X.matrix() - result.Ec.X.matrix() ).norm(), 0.0001);
    EXPECT_NEAR(result.eta, 1.2450, 1E-4);
    EXPECT_NEAR(result.omega, 3.0920, 1E-4);


}
int main(int argc, char* argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

