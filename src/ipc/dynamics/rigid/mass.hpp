#pragma once

#include <ipc/utils/eigen_ext.hpp>

#include <Eigen/Core>
#include <Eigen/Sparse>

namespace ipc::rigid {

/// @brief Compute the total mass, center of mass, and moment of inertia
/// @param vertices Vertices of the mesh
/// @param facets Facets (2D: edges; 3D: triangles) of the mesh
/// @param total_mass Total mass of the mesh
/// @param center_of_mass Center of mass of the mesh
/// @param moment_of_inertia Moment of inertia of the mesh
void compute_mass_properties(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& facets,
    double& total_mass,
    VectorMax3d& center_of_mass,
    MatrixMax3d& moment_of_inertia);

/// @brief Construct the sparse mass matrix for the given mesh (V, E).
/// @param vertices Vertices of the mesh
/// @param facets Facets (2D: edges; 3D: triangles) of the mesh
/// @param mass_matrix Sparse mass matrix of the mesh
void construct_mass_matrix(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& facets,
    Eigen::SparseMatrix<double>& mass_matrix);

} // namespace ipc::rigid
