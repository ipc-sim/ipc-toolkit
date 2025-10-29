#pragma once

#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

/// @brief Compute the bending angle between two triangles sharing an edge.
///     x0---x2
///      | \ |
///     x1---x3
/// @param x0 The first vertex of the edge.
/// @param x1 The second vertex of the edge.
/// @param x2 The opposite vertex of the first triangle.
/// @param x3 The opposite vertex of the second triangle.
/// @return The bending angle between the two triangles.
double dihedral_angle(
    Eigen::ConstRef<Eigen::Vector3d> x0,
    Eigen::ConstRef<Eigen::Vector3d> x1,
    Eigen::ConstRef<Eigen::Vector3d> x2,
    Eigen::ConstRef<Eigen::Vector3d> x3);

/// @brief Compute the Jacobian of the bending angle between two triangles sharing an edge.
///     x0---x2
///      | \ |
///     x1---x3
/// @param x0 The first vertex of the edge.
/// @param x1 The second vertex of the edge.
/// @param x2 The opposite vertex of the first triangle.
/// @param x3 The opposite vertex of the second triangle.
/// @return The Jacobian matrix of the bending angle with respect to the input vertices.
Eigen::Vector<double, 12> dihedral_angle_gradient(
    Eigen::ConstRef<Eigen::Vector3d> x0,
    Eigen::ConstRef<Eigen::Vector3d> x1,
    Eigen::ConstRef<Eigen::Vector3d> x2,
    Eigen::ConstRef<Eigen::Vector3d> x3);

} // namespace ipc