#pragma once

#include <Eigen/Dense>
#include <ceres/ceres.h>
#include <cassert>
#include <array>
namespace scd{
    template<uint n>
    struct Superquadric{
        Superquadric() = default;
        std::array<double,n> a;
        std::array<double,n-1> e;
        Eigen::Transform<double,n,Eigen::Isometry> X;
    };

    template<uint n>
    struct CollideRequest{
        CollideRequest(){
            // Default Ceres solver options for scd
            ceres_options.minimizer_type = ceres::TRUST_REGION;
            ceres_options.max_num_iterations = 100;
            ceres_options.linear_solver_type = ceres::DENSE_QR;
            ceres_options.minimizer_progress_to_stdout = false;
        }
        ceres::Solver::Options ceres_options;
    };

    template<uint n>
    struct CollideResult{
        bool collide;
        std::array<double,n-1> angles;
        Superquadric<n> Ec;
        Eigen::Matrix<double,n,1> x_eb;
        ceres::Solver::Summary ceres_summary;
    };

    struct Request{
        Request(){
            // Default Ceres solver options for scd
            ceres_options.minimizer_type = ceres::TRUST_REGION;
            ceres_options.max_num_iterations = 100;
            ceres_options.linear_solver_type = ceres::DENSE_QR;
            ceres_options.minimizer_progress_to_stdout = false;
        }
        ceres::Solver::Options ceres_options;
    };
    struct Request2{
        Request2(){
            // Default Ceres solver options for scd
            ceres_options.minimizer_type = ceres::TRUST_REGION;
            ceres_options.max_num_iterations = 100;
            ceres_options.linear_solver_type = ceres::DENSE_QR;
            ceres_options.minimizer_progress_to_stdout = false;
        }
        ceres::Solver::Options ceres_options;
    };

    void Collide(const Superquadric<2>& SQ1, const Superquadric<2>& E2, const CollideRequest<2>& request, CollideResult<2>& result);
    void Collide(const Superquadric<3>& SQ1, const Superquadric<3>& E2, const CollideRequest<3>& request, CollideResult<3>& result);

} //namespace scd
