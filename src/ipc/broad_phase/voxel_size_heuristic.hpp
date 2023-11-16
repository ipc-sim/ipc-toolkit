#pragma once

#include <Eigen/Core>

namespace ipc {

double suggest_good_voxel_size(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& edges,
    const double inflation_radius = 0);

double suggest_good_voxel_size(
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const Eigen::MatrixXi& edges,
    const double inflation_radius = 0);

/// @brief Compute the average edge length of a mesh.
double mean_edge_length(
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const Eigen::MatrixXi& edges,
    double& std_deviation);

/// @brief Compute the average displacement length.
double mean_displacement_length(
    const Eigen::MatrixXd& displacements, double& std_deviation);

/// @brief Compute the median edge length of a mesh.
double median_edge_length(
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const Eigen::MatrixXi& edges);

/// @brief Compute the median displacement length.
double median_displacement_length(const Eigen::MatrixXd& displacements);

/// @brief Compute the maximum edge length of a mesh.
double max_edge_length(
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const Eigen::MatrixXi& edges);

/// @brief Compute the maximum displacement length.
double max_displacement_length(const Eigen::MatrixXd& displacements);

} // namespace ipc
