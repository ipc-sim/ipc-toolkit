#pragma once

#include <catch2/catch_test_macros.hpp>

#include <ipc/collisions/tangential/tangential_collision.hpp>
#include <ipc/smooth_contact/smooth_collisions.hpp>

struct FrictionData {
    Eigen::MatrixXd V0;
    Eigen::MatrixXd V1;
    Eigen::MatrixXi E;
    Eigen::MatrixXi F;
    ipc::NormalCollisions collisions;
    double mu;
    double epsv_times_h;
    double p;
    double barrier_stiffness;
};

Eigen::VectorXd LogSpaced(int num, double start, double stop, double base = 10);
Eigen::VectorXd GeomSpaced(int num, double start, double stop);

FrictionData friction_data_generator();

struct SmoothFrictionData {
    Eigen::MatrixXd V0;
    Eigen::MatrixXd V1;
    Eigen::MatrixXi E;
    Eigen::MatrixXi F;
    ipc::SmoothCollisions collisions;
    double mu;
    double epsv_times_h;
    ipc::SmoothContactParameters p;
    double barrier_stiffness;
};

SmoothFrictionData smooth_friction_data_generator_2d();
SmoothFrictionData smooth_friction_data_generator_3d();
