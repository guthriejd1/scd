#include <gtest/gtest.h>
#include <cmath>
#include <fstream>
#include "../include/scd/scd.h"
#include "data_utilities.h"
#include <iostream>

GTEST_TEST(Collision, BoxEllipse2D){
    Eigen::Rotation2D<double> Xr(-M_PI/5);
    Eigen::Translation<double,2> Xt(1,-1);

    scd::Superquadric<2> SQ1;
    SQ1.X = Xt*Xr;
    SQ1.a = {1.0,2.0};
    SQ1.e = {0.25};

    Xr.angle() = -M_PI/3;
    Xt.x() = -3;
    Xt.y() = 3;

    scd::Superquadric<2> E2;
    E2.X = Xt*Xr;
    E2.a = {0.5,2.0};
    E2.e = {1};

    scd::CollideRequest<2> request;
    scd::CollideResult<2> result;
    scd::Collide(SQ1,E2,request,result);

    Eigen::Transform<double,2,Eigen::Isometry> X_expect;
    X_expect.matrix() <<  0.5000,    0.8660,   -0.3587,
                         -0.8660,    0.5000,    0.3587,
                          0,         0,         1.0000;
    ASSERT_LT(( X_expect.matrix() - result.Ec.X.matrix() ).norm(), 1E-4);
    EXPECT_NEAR(result.angles[0],2.726233503871522, 1E-6);
    ASSERT_FALSE(result.collide);
}

GTEST_TEST(Collision, CubeEllipsoid3D){
    auto Xr = Eigen::AngleAxis<double>(M_PI/3, Eigen::Vector3d::UnitX())
                 * Eigen::AngleAxis<double>(2*M_PI/3, Eigen::Vector3d::UnitY());
    auto Xt = Eigen::Translation<double,3>(-1.0,1.0,2.0);

    scd::Superquadric<3> SQ1;
    SQ1.X = Xt*Xr;
    SQ1.a = {1.0,2.0,3.0};
    SQ1.e = {0.25,0.25};

    Xr = Eigen::AngleAxis<double>(7*M_PI/6, Eigen::Vector3d::UnitX())
                  * Eigen::AngleAxis<double>(-4*M_PI/6, Eigen::Vector3d::UnitY())
                  * Eigen::AngleAxis<double>(3*M_PI/5, Eigen::Vector3d::UnitZ());
    Xt = Eigen::Translation<double,3>(9.0,3.0,6.0);

    scd::Superquadric<3> E2;
    E2.X = Xt*Xr;
    E2.a = {1.0,2.0,0.5};
    E2.e = {1.0,1.0};

    scd::CollideRequest<3> request;
    scd::CollideResult<3> result;
    scd::Collide(SQ1,E2,request,result);

    Eigen::Transform<double,3,Eigen::Isometry> X;
    X.matrix() <<      0.1545,    0.4755,   -0.8660,    2.9909,
                      -0.9574,   -0.1442,   -0.2500,    1.7982,
                      -0.2438,    0.8678,    0.4330,    3.5964,
                       0,         0,         0,         1.0000;
    // Data was copied from MATLAB printout, keeping only 4 decimals. Tolerances reflect this.

    ASSERT_LT(( X.matrix() - result.Ec.X.matrix() ).norm(), 0.0001);
    EXPECT_NEAR(result.angles[0], 3.0920, 1E-4);
    EXPECT_NEAR(result.angles[1], 1.2450, 1E-4);
}

GTEST_TEST(Collision, EllipsoidEllipsoidBatch){
    std::string file_path = "../test/data/ellipsoid_ellipsoid_batch.txt";
    uint n_line = 0;
    std::cout << "Begin batch test...\n";
    for (int n_test = 0; n_test < 1000; n_test++) {
        scd::Superquadric<3> SQ1;
        scd::Superquadric<3> E2;
        scd::Superquadric<3> E2c_expect;
        scd::CollideResult<3> result_expect;
        scd::test::ReadTestFile(file_path, n_line, SQ1, E2, E2c_expect, result_expect);
        scd::CollideRequest<3> request;
        scd::CollideResult<3> result;

        scd::Collide(SQ1,E2,request,result);

       ASSERT_EQ(result.collide, result_expect.collide);
       // TODO (jay) Add angle wrapping to allow comparison
        n_line += 25;
    }
}

int main(int argc, char* argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

