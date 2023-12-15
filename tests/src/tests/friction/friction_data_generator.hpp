#pragma once

#include <catch2/catch_test_macros.hpp>

#include <ipc/friction/collisions/friction_collision.hpp>

struct FrictionData {
    Eigen::MatrixXd V0;
    Eigen::MatrixXd V1;
    Eigen::MatrixXi E;
    Eigen::MatrixXi F;
    ipc::Collisions collisions;
    double mu;
    double epsv_times_h;
    double p;
    double barrier_stiffness;
};

Eigen::VectorXd LogSpaced(int num, double start, double stop, double base = 10);
Eigen::VectorXd GeomSpaced(int num, double start, double stop);

FrictionData friction_data_generator();