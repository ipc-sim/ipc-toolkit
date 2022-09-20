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
double mean_edge_length(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    double& std_deviation);

/// @brief Compute the average displacement length.
double
mean_displacement_length(const Eigen::MatrixXd& U, double& std_deviation);

/// @brief Compute the median edge length of a mesh.
double median_edge_length(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E);

/// @brief Compute the median displacement length.
double median_displacement_length(const Eigen::MatrixXd& U);

} // namespace ipc
