#pragma once

#include <catch2/catch_test_macros.hpp>

#include <ipc/collisions/tangential/tangential_collision.hpp>
#include <ipc/gcp/gcp_collisions.hpp>

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

struct GCPFrictionData {
    Eigen::MatrixXd V0;
    Eigen::MatrixXd V1;
    Eigen::MatrixXi E;
    Eigen::MatrixXi F;
    ipc::GCPCollisions collisions;
    double mu;
    double epsv_times_h;
    ipc::GCPParameters p;
    double barrier_stiffness;
};

GCPFrictionData smooth_friction_data_generator_2d();
GCPFrictionData smooth_friction_data_generator_3d();

/// Scene geometry for "ESP friction force jacobian 3D" tests.
/// Sections: "point-triangle", "point-edge", "point-point".
struct ESPFrictionSceneData3D {
    Eigen::MatrixXd X;
    Eigen::MatrixXi E;
    Eigen::MatrixXi F;
    /// Vertex indices of the "upper object" that slides during the test.
    std::vector<int> upper_vertices;
};

ESPFrictionSceneData3D esp_friction_scene_generator_3d(double d);
