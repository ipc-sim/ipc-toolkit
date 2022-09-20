#pragma once

#include <catch2/catch.hpp>

#include <Eigen/Core>

#include <ipc/friction/friction_constraint.hpp>

struct FrictionData {
    Eigen::MatrixXd V0;
    Eigen::MatrixXd V1;
    Eigen::MatrixXi E;
    Eigen::MatrixXi F;
    ipc::Constraints constraints;
    double mu;
    double epsv_times_h;
    double dhat;
    double barrier_stiffness;
};

Eigen::VectorXd LogSpaced(int num, double start, double stop, double base = 10);
Eigen::VectorXd GeomSpaced(int num, double start, double stop);

FrictionData friction_data_generator();
