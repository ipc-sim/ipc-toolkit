#pragma once

#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

/// @brief Suggest a good voxel size for the given mesh.
/// @param vertices The vertex positions of the mesh.
/// @param edges The edges of the mesh.
/// @param inflation_radius The radius of inflation around all elements.
/// @return The suggested voxel size.
double suggest_good_voxel_size(
    Eigen::ConstRef<Eigen::MatrixXd> vertices,
    Eigen::ConstRef<Eigen::MatrixXi> edges,
    const double inflation_radius = 0);

/// @brief Suggest a good voxel size for the given mesh.
/// @param vertices_t0 The vertex positions of the mesh at time t0.
/// @param vertices_t1 The vertex positions of the mesh at time t1.
/// @param edges The edges of the mesh.
/// @param inflation_radius The radius of inflation around all elements.
/// @return The suggested voxel size.
double suggest_good_voxel_size(
    Eigen::ConstRef<Eigen::MatrixXd> vertices_t0,
    Eigen::ConstRef<Eigen::MatrixXd> vertices_t1,
    Eigen::ConstRef<Eigen::MatrixXi> edges,
    const double inflation_radius = 0);

/// @brief Compute the average edge length of a mesh.
/// @param vertices_t0 The vertex positions of the mesh at time t0.
/// @param vertices_t1 The vertex positions of the mesh at time t1.
/// @param edges The edges of the mesh.
/// @param[out] std_deviation The standard deviation of the edge lengths.
/// @return The average edge length.
double mean_edge_length(
    Eigen::ConstRef<Eigen::MatrixXd> vertices_t0,
    Eigen::ConstRef<Eigen::MatrixXd> vertices_t1,
    Eigen::ConstRef<Eigen::MatrixXi> edges,
    double& std_deviation);

/// @brief Compute the average displacement length.
/// @param displacements The displacements of the vertices.
/// @param[out] std_deviation The standard deviation of the displacement lengths.
/// @return The average displacement length.
double mean_displacement_length(
    Eigen::ConstRef<Eigen::MatrixXd> displacements, double& std_deviation);

/// @brief Compute the median edge length of a mesh.
/// @param vertices_t0 The vertex positions of the mesh at time t0.
/// @param vertices_t1 The vertex positions of the mesh at time t1.
/// @param edges The edges of the mesh.
/// @return The median edge length.
double median_edge_length(
    Eigen::ConstRef<Eigen::MatrixXd> vertices_t0,
    Eigen::ConstRef<Eigen::MatrixXd> vertices_t1,
    Eigen::ConstRef<Eigen::MatrixXi> edges);

/// @brief Compute the median displacement length.
/// @param displacements The displacements of the vertices.
double
median_displacement_length(Eigen::ConstRef<Eigen::MatrixXd> displacements);

/// @brief Compute the maximum edge length of a mesh.
/// @param vertices_t0 The vertex positions of the mesh at time t0.
/// @param vertices_t1 The vertex positions of the mesh at time t1.
/// @param edges The edges of the mesh.
/// @return The maximum edge length.
double max_edge_length(
    Eigen::ConstRef<Eigen::MatrixXd> vertices_t0,
    Eigen::ConstRef<Eigen::MatrixXd> vertices_t1,
    Eigen::ConstRef<Eigen::MatrixXi> edges);

/// @brief Compute the maximum displacement length.
/// @param displacements The displacements of the vertices.
/// @return The maximum displacement length.
double max_displacement_length(Eigen::ConstRef<Eigen::MatrixXd> displacements);

} // namespace ipc
