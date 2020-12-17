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
            Superquadric3 SQ1, Ellipsoid3 E2){
        // Decision variables
        Eigen::Map<const Eigen::Matrix<T,1,1>> omega(omega_ptr);
        Eigen::Map<const Eigen::Matrix<T,1,1>> eta(eta_ptr);

        double r = std::min(std::min(E2.ax,E2.ay),E2.az);
        Eigen::Matrix<double,3,3> E2_R = E2.X.matrix().block<3,3>(0,0);
        Eigen::Matrix<double,3,1> E2_t = E2.X.matrix().block<3,1>(0,3);

        Eigen::Matrix<double,3,1> d;
        d << r/E2.ax, r/E2.ay, r/E2.az;
        Eigen::Matrix<double,3,3> M = E2_R*d.asDiagonal()*E2_R.transpose();

        double omega_cos_sign = 1.0;
        double omega_sin_sign = 1.0;
        auto omega_mod = ModifyAngleQuadrant(omega[0], omega_cos_sign, omega_sin_sign);

        double eta_cos_sign = 1.0;
        double eta_sin_sign = 1.0;
        auto eta_mod = ModifyAngleQuadrant(eta[0], eta_cos_sign, eta_sin_sign);

        Eigen::Matrix<T,3,1> x_a;
        x_a(0) = SQ1.ax*eta_cos_sign*omega_cos_sign*pow(cos(eta_mod),SQ1.e1)*pow(cos(omega_mod),SQ1.e2);
        x_a(1) = SQ1.ay*eta_cos_sign*omega_sin_sign*pow(cos(eta_mod),SQ1.e1)*pow(sin(omega_mod),SQ1.e2);
        x_a(2) = SQ1.az*eta_sin_sign*pow(sin(eta_mod),SQ1.e1);

        Eigen::Matrix<T,3,1> grad_phi;
        grad_phi(0) = 1/SQ1.ax*eta_cos_sign*omega_cos_sign*pow(cos(eta_mod),2.0-SQ1.e1)*pow(cos(omega_mod),2.0-SQ1.e2);
        grad_phi(1) = 1/SQ1.ay*eta_cos_sign*omega_sin_sign*pow(cos(eta_mod),2.0-SQ1.e1)*pow(sin(omega_mod),2.0-SQ1.e2);
        grad_phi(2) = 1/SQ1.az*eta_sin_sign*pow(sin(eta_mod),2.0-SQ1.e1);

        Eigen::Matrix<T,3,1> x_eb;
        Eigen::Matrix<T,3,1> denom = M.inverse()*grad_phi;
        x_eb = x_a + r*(M.inverse()*M.inverse())*grad_phi / denom.norm() ;
        return x_eb;
    }
    template <typename T>
    Eigen::Matrix<T,2,1> Calculate_x_eb2(const T* const theta_ptr,
                                        Superquadric2 SQ1, Superquadric2 E2){
        // Decision variables
        Eigen::Map<const Eigen::Matrix<T,1,1>> theta(theta_ptr);

        double r = std::min(E2.ax,E2.ay);
        Eigen::Matrix<double,2,2> E2_R = E2.X.matrix().block<2,2>(0,0);
        Eigen::Matrix<double,2,1> E2_t = E2.X.matrix().block<2,1>(0,2);

        Eigen::Matrix<double,2,1> d;
        d << r/E2.ax, r/E2.ay;
        Eigen::Matrix<double,2,2> M = E2_R*d.asDiagonal()*E2_R.transpose();

        double theta_cos_sign = 1.0;
        double theta_sin_sign = 1.0;
        auto theta_mod = ModifyAngleQuadrant(theta[0], theta_cos_sign, theta_sin_sign);

        Eigen::Matrix<T,2,1> x_a;
        x_a(0) = SQ1.ax*theta_cos_sign*pow(cos(theta_mod),SQ1.e);
        x_a(1) = SQ1.ay*theta_sin_sign*pow(sin(theta_mod),SQ1.e);

        Eigen::Matrix<T,2,1> grad_phi;
        grad_phi(0) = 2/SQ1.ax*theta_cos_sign*pow(cos(theta_mod),2.0-SQ1.e);
        grad_phi(1) = 2/SQ1.ay*theta_sin_sign*pow(sin(theta_mod),2.0-SQ1.e);

        Eigen::Matrix<T,2,1> x_eb;
        Eigen::Matrix<T,2,1> denom = M.inverse()*grad_phi;
        x_eb = x_a + r*(M.inverse()*M.inverse())*grad_phi / denom.norm() ;
        return x_eb;
    }

    struct RootResidual3 {
        RootResidual3(Superquadric3 m_SQ1, Ellipsoid3 m_E2)
                : SQ1(m_SQ1), E2(m_E2){}
        template <typename T>
        bool operator()(const T* const omega_ptr, const T* const eta_ptr, T* res_ptr) const {
            // Residuals
            Eigen::Map<Eigen::Matrix<T,3,1>> res(res_ptr);

            Eigen::Matrix<T,3,1> x_eb = Calculate_x_eb(omega_ptr, eta_ptr, SQ1, E2);
            res = x_eb.cross(E2.X.matrix().block<3,1>(0,3));

            return true;
        }
    private:
        Superquadric3 SQ1;
        Ellipsoid3 E2;
    };

    struct RootResidual2 {
        RootResidual2(Superquadric2 m_SQ1, Superquadric2 m_E2)
                : SQ1(m_SQ1), E2(m_E2){}
        template <typename T>
        bool operator()(const T* const theta_ptr, T* res_ptr) const {
            // Residuals
            Eigen::Map<Eigen::Matrix<T,1,1>> res(res_ptr);

            Eigen::Matrix<T,2,1> x_eb = Calculate_x_eb2(theta_ptr, SQ1, E2);
            res = E2.X.matrix().block<1,1>(1,2)*x_eb(0) - E2.X.matrix().block<1,1>(0,2)*x_eb(1);

            return true;
        }
    private:
        Superquadric2 SQ1;
        Superquadric2 E2;
    };

    void collide(const Superquadric3& SQ1_W, const Ellipsoid3& E2_W, const Request& request, Result& result){
        // Transform SQ1 and E2 from world frame to SQ1 frame
        Superquadric3 SQ1 = SQ1_W;
        SQ1.X.setIdentity();
        Ellipsoid3 E2 = E2_W;
        E2.X = SQ1_W.X.inverse()*E2_W.X;

        // TODO(jay) Improve the initial guess code
        // Generate initial guess for decision variables omega, eta
        double omega = atan2(E2.X(1,3),E2.X(0,3));
        if (omega < 0) {omega = omega + 2*M_PI;}

        double eta = 0.0;
        if (E2.X(2,3) >= 0)
            eta = M_PI/4;
        else
            eta = 2.0*M_PI - M_PI/4;

        // Setup Ceres problem
        ceres::Problem problem;
        problem.AddResidualBlock(
                new ceres::AutoDiffCostFunction<RootResidual3, 3, 1, 1>(new RootResidual3(SQ1, E2)), NULL, &omega, &eta);

        ceres::Solve(request.ceres_options, &problem, &result.ceres_summary);

        // TODO(jay) Add option for multiple solve attempts

        result.omega = omega;
        result.eta = eta;

        // Determine if superquadric and ellipsoid are in collision
        result.x_eb = Calculate_x_eb(&result.omega, &result.eta, SQ1, E2);
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

    void collide2(const Superquadric2& SQ1_W, const Superquadric2& E2_W, const Request2& request, Result2& result){
        // Transform SQ1 and E2 from world frame to SQ1 frame
        Superquadric2 SQ1 = SQ1_W;
        SQ1.X.setIdentity();
        Superquadric2 E2 = E2_W;
        E2.X = SQ1_W.X.inverse()*E2_W.X;

        // Generate initial guess for decision variable theta
        double theta = atan2(E2.X(1,2),E2.X(0,2));
        if (theta < 0) {theta = theta + 2*M_PI;}

        // Setup Ceres problem
        ceres::Problem problem;
        problem.AddResidualBlock(
                new ceres::AutoDiffCostFunction<RootResidual2, 1, 1>(new RootResidual2(SQ1, E2)), NULL, &theta);

        ceres::Solve(request.ceres_options, &problem, &result.ceres_summary);

        result.theta = theta;

        // Determine if superquadric and ellipsoid are in collision
        result.x_eb = Calculate_x_eb2(&result.theta, SQ1, E2);
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
}// namespace scd

