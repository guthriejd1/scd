#include <gtest/gtest.h>
#include <cmath>>
#include "../include/scd/scd.h"


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

    // TODO(jay) Load validation data from file to allow more numerical precision
    Eigen::Transform<double,3,Eigen::Isometry> X;
    X.matrix() <<      0.1545,    0.4755,   -0.8660,    2.9909,
                      -0.9574,   -0.1442,   -0.2500,    1.7982,
                      -0.2438,    0.8678,    0.4330,    3.5964,
                       0,         0,         0,         1.0000;

    ASSERT_LT(( X.matrix() - result.Ec.X.matrix() ).norm(), 0.0001);
    EXPECT_NEAR(result.eta, 1.2450, 1E-4);
    EXPECT_NEAR(result.omega, 3.0920, 1E-4);
}
int main(int argc, char* argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

