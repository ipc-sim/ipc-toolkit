#pragma once

#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

// ============================================================================
// Point - Point

/// @brief Compute the relative velocity of two points
/// @param dp0 Velocity of the first point
/// @param dp1 Velocity of the second point
/// @return The relative velocity of the two points
VectorMax3d point_point_relative_velocity(
    const Eigen::Ref<const VectorMax3d>& dp0,
    const Eigen::Ref<const VectorMax3d>& dp1);

/// @brief Compute the relative velocity premultiplier matrix
/// @param dim Dimension (2 or 3)
/// @return The relative velocity premultiplier matrix
MatrixMax<double, 3, 6> point_point_relative_velocity_matrix(const int dim);

/// @brief Compute the jacobian of the relative velocity premultiplier matrix
/// @param dim Dimension (2 or 3)
/// @return The jacobian of the relative velocity premultiplier matrix
MatrixMax<double, 3, 6>
point_point_relative_velocity_matrix_jacobian(const int dim);

// ============================================================================
// Point - Edge

/// @brief Compute the relative velocity of a point and an edge
/// @param dp Velocity of the point
/// @param de0 Velocity of the first endpoint of the edge
/// @param de1 Velocity of the second endpoint of the edge
/// @param alpha Parametric coordinate of the closest point on the edge
/// @return The relative velocity of the point and the edge
VectorMax3d point_edge_relative_velocity(
    const Eigen::Ref<const VectorMax3d>& dp,
    const Eigen::Ref<const VectorMax3d>& de0,
    const Eigen::Ref<const VectorMax3d>& de1,
    const double alpha);

MatrixMax<double, 3, 9>
point_edge_relative_velocity_matrix(const int dim, const double alpha);

MatrixMax<double, 3, 9>
point_edge_relative_velocity_matrix_jacobian(const int dim, const double alpha);

// ============================================================================
// Edge - Edge

/// @brief Compute the relative velocity of the edges.
/// @param dea0 Velocity of the first endpoint of the first edge
/// @param dea1 Velocity of the second endpoint of the first edge
/// @param deb0 Velocity of the first endpoint of the second edge
/// @param deb1 Velocity of the second endpoint of the second edge
/// @param coords Two parametric coordinates of the closest points on the edges
/// @return The relative velocity of the edges
Eigen::Vector3d edge_edge_relative_velocity(
    const Eigen::Ref<const Eigen::Vector3d>& dea0,
    const Eigen::Ref<const Eigen::Vector3d>& dea1,
    const Eigen::Ref<const Eigen::Vector3d>& deb0,
    const Eigen::Ref<const Eigen::Vector3d>& deb1,
    const Eigen::Ref<const Eigen::Vector2d>& coords);

MatrixMax<double, 3, 12> edge_edge_relative_velocity_matrix(
    const int dim, const Eigen::Ref<const Eigen::Vector2d>& coords);

MatrixMax<double, 6, 12> edge_edge_relative_velocity_matrix_jacobian(
    const int dim, const Eigen::Ref<const Eigen::Vector2d>& coords);

// ============================================================================
// Point - Triangle

/// @brief Compute the relative velocity of the point to the triangle.
/// @param dp Velocity of the point
/// @param dt0 Velocity of the first vertex of the triangle
/// @param dt1 Velocity of the second vertex of the triangle
/// @param dt2 Velocity of the third vertex of the triangle
/// @param coords Baricentric coordinates of the closest point on the triangle
/// @return The relative velocity of the point to the triangle
Eigen::Vector3d point_triangle_relative_velocity(
    const Eigen::Ref<const Eigen::Vector3d>& dp,
    const Eigen::Ref<const Eigen::Vector3d>& dt0,
    const Eigen::Ref<const Eigen::Vector3d>& dt1,
    const Eigen::Ref<const Eigen::Vector3d>& dt2,
    const Eigen::Ref<const Eigen::Vector2d>& coords);

MatrixMax<double, 3, 12> point_triangle_relative_velocity_matrix(
    const int dim, const Eigen::Ref<const Eigen::Vector2d>& coords);

MatrixMax<double, 6, 12> point_triangle_relative_velocity_matrix_jacobian(
    const int dim, const Eigen::Ref<const Eigen::Vector2d>& coords);

} // namespace ipc
