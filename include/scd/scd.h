#pragma once

#include <Eigen/Dense>
#include <ceres/ceres.h>
#include <cassert>

namespace scd{
    struct Ellipsoid3 {
        Ellipsoid3(double m_ax, double m_ay, double m_az) : ax(m_ax), ay(m_ay), az(m_az) {
            assert((m_ax > 0) && (m_ay > 0) && (m_az > 0));
        }
        Ellipsoid3() = default;
        double ax;
        double ay;
        double az;
        Eigen::Transform<double,3,Eigen::Isometry> X;
    };

    struct Superquadric3{
        Superquadric3(double m_ax, double m_ay, double m_az, double m_e1, double m_e2) : ax(m_ax), ay(m_ay), az(m_az), e1(m_e1), e2(m_e2) {
            assert( (m_ax > 0) && (m_ay > 0) && (m_az > 0) );
            assert( (m_e1 >= 0) && (m_e1 <= 2) );
            assert( (m_e2 >= 0) && (m_e2 <= 2) );
        }
        double ax;
        double ay;
        double az;
        double e1;
        double e2;
        Eigen::Transform<double,3,Eigen::Isometry> X;
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

    struct Result{
        bool collide;
        double omega;
        double eta;
        Ellipsoid3 Ec;
        Eigen::Matrix<double,3,1> x_eb;
        ceres::Solver::Summary ceres_summary;
    };

    void collide(const Superquadric3& SQ1, const Ellipsoid3& E2, const Request& request, Result& result);

} //namespace scd
