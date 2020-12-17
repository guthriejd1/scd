#include <iostream>

#include <Eigen/Dense>
#include <ceres/ceres.h>

#include "../include/scd/scd.h"

namespace scd{


    template <typename T>
    T ModifyAngleQuadrant(const T& theta, double& cos_sign, double& sin_sign){
        // TODO(jay) Add support for angles outside [0 2pi] interval
        auto theta_mod = theta;
        if ( (theta >= 0.0) && (theta <= M_PI/2) ){
        }
        else if ( (theta >= M_PI/2) && (theta <= M_PI) ){
            theta_mod = M_PI - theta;
            cos_sign = -1.0;
        }
        else if ( (theta >= M_PI) && (theta <= 1.5*M_PI) ){
            theta_mod = -M_PI + theta;
            cos_sign = -1.0;
            sin_sign = -1.0;
        }
        else if ( (theta >= 1.5*M_PI) && (theta <= 2.0*M_PI) ){
            theta_mod = 2.0*M_PI - theta;
            sin_sign = -1.0;
        }
        return theta_mod;
    }

    template <typename T>
    Eigen::Matrix<T,3,1> Calculate_x_eb(const T* const omega_ptr, const T* const eta_ptr,
            Superquadric<3> SQ1, Superquadric<3> E2){
        // Decision variables
        Eigen::Map<const Eigen::Matrix<T,1,1>> omega(omega_ptr);
        Eigen::Map<const Eigen::Matrix<T,1,1>> eta(eta_ptr);

        double r = std::min(std::min(E2.a[0],E2.a[1]),E2.a[2]);
        Eigen::Matrix<double,3,3> E2_R = E2.X.matrix().block<3,3>(0,0);
        Eigen::Matrix<double,3,1> E2_t = E2.X.matrix().block<3,1>(0,3);

        Eigen::Matrix<double,3,1> d;
        d << r/E2.a[0], r/E2.a[1], r/E2.a[2];
        Eigen::Matrix<double,3,3> M = E2_R*d.asDiagonal()*E2_R.transpose();
        double omega_cos_sign = 1.0;
        double omega_sin_sign = 1.0;
        auto omega_mod = ModifyAngleQuadrant(omega[0], omega_cos_sign, omega_sin_sign);

        double eta_cos_sign = 1.0;
        double eta_sin_sign = 1.0;
        auto eta_mod = ModifyAngleQuadrant(eta[0], eta_cos_sign, eta_sin_sign);

        Eigen::Matrix<T,3,1> x_a;
        x_a(0) = SQ1.a[0]*eta_cos_sign*omega_cos_sign*pow(cos(eta_mod),SQ1.e[0])*pow(cos(omega_mod),SQ1.e[1]);
        x_a(1) = SQ1.a[1]*eta_cos_sign*omega_sin_sign*pow(cos(eta_mod),SQ1.e[0])*pow(sin(omega_mod),SQ1.e[1]);
        x_a(2) = SQ1.a[2]*eta_sin_sign*pow(sin(eta_mod),SQ1.e[0]);

        Eigen::Matrix<T,3,1> grad_phi;
        grad_phi(0) = 1/SQ1.a[0]*eta_cos_sign*omega_cos_sign*pow(cos(eta_mod),2.0-SQ1.e[0])*pow(cos(omega_mod),2.0-SQ1.e[1]);
        grad_phi(1) = 1/SQ1.a[1]*eta_cos_sign*omega_sin_sign*pow(cos(eta_mod),2.0-SQ1.e[0])*pow(sin(omega_mod),2.0-SQ1.e[1]);
        grad_phi(2) = 1/SQ1.a[2]*eta_sin_sign*pow(sin(eta_mod),2.0-SQ1.e[0]);

        Eigen::Matrix<T,3,1> x_eb;
        Eigen::Matrix<T,3,1> denom = M.inverse()*grad_phi;
        x_eb = x_a + r*(M.inverse()*M.inverse())*grad_phi / denom.norm() ;
        return x_eb;
    }
    template <typename T>
    Eigen::Matrix<T,2,1> Calculate_x_eb(const T* const theta_ptr,
                                        Superquadric<2> SQ1, Superquadric<2> E2){
        // Decision variables
        Eigen::Map<const Eigen::Matrix<T,1,1>> theta(theta_ptr);

        double r = std::min(E2.a[0],E2.a[1]);
        Eigen::Matrix<double,2,2> E2_R = E2.X.matrix().block<2,2>(0,0);
        Eigen::Matrix<double,2,1> E2_t = E2.X.matrix().block<2,1>(0,2);

        Eigen::Matrix<double,2,1> d;
        d << r/E2.a[0], r/E2.a[1];
        Eigen::Matrix<double,2,2> M = E2_R*d.asDiagonal()*E2_R.transpose();

        double theta_cos_sign = 1.0;
        double theta_sin_sign = 1.0;
        auto theta_mod = ModifyAngleQuadrant(theta[0], theta_cos_sign, theta_sin_sign);

        Eigen::Matrix<T,2,1> x_a;
        x_a(0) = SQ1.a[0]*theta_cos_sign*pow(cos(theta_mod),SQ1.e[0]);
        x_a(1) = SQ1.a[1]*theta_sin_sign*pow(sin(theta_mod),SQ1.e[0]);

        Eigen::Matrix<T,2,1> grad_phi;
        grad_phi(0) = 2/SQ1.a[0]*theta_cos_sign*pow(cos(theta_mod),2.0-SQ1.e[0]);
        grad_phi(1) = 2/SQ1.a[1]*theta_sin_sign*pow(sin(theta_mod),2.0-SQ1.e[0]);

        Eigen::Matrix<T,2,1> x_eb;
        Eigen::Matrix<T,2,1> denom = M.inverse()*grad_phi;
        x_eb = x_a + r*(M.inverse()*M.inverse())*grad_phi / denom.norm() ;
        return x_eb;
    }


    struct Residual2 {
        Residual2(Superquadric<2> m_SQ1, Superquadric<2> m_E2)
                : SQ1(m_SQ1), E2(m_E2){}
        template <typename T>
        bool operator()(const T* const theta_ptr, T* res_ptr) const {
            // Residuals
            Eigen::Map<Eigen::Matrix<T,1,1>> res(res_ptr);

            Eigen::Matrix<T,2,1> x_eb = Calculate_x_eb(theta_ptr, SQ1, E2);
            res = E2.X.matrix().block<1,1>(1,2)*x_eb(0) - E2.X.matrix().block<1,1>(0,2)*x_eb(1);
            return true;
        }
    private:
        Superquadric<2> SQ1;
        Superquadric<2> E2;
    };


    struct Residual3 {
        Residual3(Superquadric<3> m_SQ1, Superquadric<3> m_E2)
                : SQ1(m_SQ1), E2(m_E2){}
        template <typename T>
        bool operator()(const T* const omega_ptr, const T* const eta_ptr, T* res_ptr) const {
            // Residuals
            Eigen::Map<Eigen::Matrix<T, 3, 1>> res(res_ptr);
            Eigen::Map<const Eigen::Matrix<T,1,1>> omega(omega_ptr);
            Eigen::Map<const Eigen::Matrix<T,1,1>> eta(eta_ptr);

            Eigen::Matrix<T, 3, 1> x_eb = Calculate_x_eb(omega_ptr, eta_ptr, SQ1, E2);
            res = x_eb.cross(E2.X.matrix().block<3, 1>(0, 3));
            return true;
        }
    private:
        Superquadric<3> SQ1;
        Superquadric<3> E2;
    };


    void Collide(const Superquadric<2>& SQ1_W, const Superquadric<2>& E2_W, const CollideRequest<2>& request, CollideResult<2>& result){
        // Transform SQ1 and E2 from world frame to SQ1 frame
        Superquadric<2> SQ1 = SQ1_W;
        SQ1.X.setIdentity();
        Superquadric<2> E2 = E2_W;
        E2.X = SQ1_W.X.inverse()*E2_W.X;

        ceres::Problem problem;
        // Generate initial guess and setup Ceres problem
        double theta = atan2(E2.X(1,2),E2.X(0,2));
        if (theta < 0) {theta = theta + 2*M_PI;}

        problem.AddResidualBlock(
                new ceres::AutoDiffCostFunction<Residual2, 1, 1>(new Residual2(SQ1, E2)), NULL, &theta);
        ceres::Solve(request.ceres_options, &problem, &result.ceres_summary);
        result.angles[0] = theta;

        // Determine if superquadric and ellipsoid are in collision
        result.x_eb = Calculate_x_eb(&result.angles[0], SQ1, E2);
        double scale = result.x_eb.norm()/E2.X.matrix().block<2,1>(0,2).norm();
        result.x_eb = scale*E2.X.matrix().block<2,1>(0,2);
        result.collide = (scale < 1.0) ? false : true;

        // Transform back to world frame to calculate pose that gives surface contact
        Eigen::Transform<double,2,Eigen::Isometry> X_SQ1_E2c;
        X_SQ1_E2c = E2.X;
        X_SQ1_E2c.matrix().block<2,1>(0,2) = result.x_eb;

        result.Ec = E2;
        result.Ec.X = SQ1_W.X*X_SQ1_E2c;
    }

    void Collide(const Superquadric<3>& SQ1_W, const Superquadric<3>& E2_W, const CollideRequest<3>& request, CollideResult<3>& result){
        // Transform SQ1 and E2 from world frame to SQ1 frame
        Superquadric<3> SQ1 = SQ1_W;
        SQ1.X.setIdentity();
        Superquadric<3> E2 = E2_W;
        E2.X = SQ1_W.X.inverse()*E2_W.X;

        ceres::Problem problem;
        // Generate initial guess and setup Ceres problem
        double omega = atan2(E2.X(1,3),E2.X(0,3));
        if (omega < 0) {omega = omega + 2*M_PI;}

        double eta = 0.0;
        if (E2.X(2,3) >= 0)
            eta = M_PI/4;
        else
            eta = 2.0*M_PI - M_PI/4;

        // Setup Ceres problem
        problem.AddResidualBlock(
                new ceres::AutoDiffCostFunction<Residual3, 3, 1, 1>(new Residual3(SQ1, E2)), NULL, &omega, &eta);

        ceres::Solve(request.ceres_options, &problem, &result.ceres_summary);

        // TODO(jay) Add option for multiple solve attempts

        result.angles[0] = omega;
        result.angles[1] = eta;

        // Determine if superquadric and ellipsoid are in collision
        result.x_eb = Calculate_x_eb(&result.angles[0], &result.angles[1], SQ1, E2);
        double scale = result.x_eb.norm()/E2.X.matrix().block<3,1>(0,3).norm();
        result.x_eb = scale*E2.X.matrix().block<3,1>(0,3);
        result.collide = (scale < 1.0) ? false : true;

        // Transform back to world frame to calculate pose that gives surface contact
        Eigen::Transform<double,3,Eigen::Isometry> X_SQ1_E2c;
        X_SQ1_E2c = E2.X;
        X_SQ1_E2c.matrix().block<3,1>(0,3) = result.x_eb;

        result.Ec = E2;
        result.Ec.X = SQ1_W.X*X_SQ1_E2c;
    }

}// namespace scd

