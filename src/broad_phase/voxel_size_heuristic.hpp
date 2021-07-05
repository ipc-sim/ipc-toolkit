#pragma once

#include <Eigen/Core>

namespace ipc {

double suggest_good_voxel_size(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    double inflation_radius = 0);

double suggest_good_voxel_size(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    double inflation_radius = 0);

/// @brief Compute the average edge length of a mesh.
double average_edge_length(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E);

/// @brief Compute the average displacement length.
double average_displacement_length(const Eigen::MatrixXd& U);

} // namespace ipc
