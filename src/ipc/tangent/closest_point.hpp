#pragma once

#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

// ============================================================================
// Point - Edge

/// @brief Compute the baricentric coordinate of the closest point on the edge.
/// @param p Point
/// @param e0 First edge point
/// @param e1 Second edge point
/// @return barycentric coordinates of the closest point
double point_edge_closest_point(
    const Eigen::Ref<const VectorMax3d>& p,
    const Eigen::Ref<const VectorMax3d>& e0,
    const Eigen::Ref<const VectorMax3d>& e1);

/// @brief Compute the Jacobian of the closest point on the edge.
/// @param p Point
/// @param e0 First edge point
/// @param e1 Second edge point
/// @return Jacobian of the closest point
VectorMax9d point_edge_closest_point_jacobian(
    const Eigen::Ref<const VectorMax3d>& p,
    const Eigen::Ref<const VectorMax3d>& e0,
    const Eigen::Ref<const VectorMax3d>& e1);

// ============================================================================
// Edge - Edge

/// @brief Compute the barycentric coordinates of the closest points between two edges.
/// @param ea0 First point of the first edge
/// @param ea1 Second point of the first edge
/// @param eb0 First point of the second edge
/// @param eb1 Second point of the second edge
/// @return Barycentric coordinates of the closest points
Eigen::Vector2d edge_edge_closest_point(
    const Eigen::Ref<const Eigen::Vector3d>& ea0,
    const Eigen::Ref<const Eigen::Vector3d>& ea1,
    const Eigen::Ref<const Eigen::Vector3d>& eb0,
    const Eigen::Ref<const Eigen::Vector3d>& eb1);

/// @brief Compute the Jacobian of the closest points between two edges.
/// @param ea0 First point of the first edge
/// @param ea1 Second point of the first edge
/// @param eb0 First point of the second edge
/// @param eb1 Second point of the second edge
/// @return Jacobian of the closest points
Eigen::Matrix<double, 2, 12> edge_edge_closest_point_jacobian(
    const Eigen::Ref<const Eigen::Vector3d>& ea0,
    const Eigen::Ref<const Eigen::Vector3d>& ea1,
    const Eigen::Ref<const Eigen::Vector3d>& eb0,
    const Eigen::Ref<const Eigen::Vector3d>& eb1);

// ============================================================================
// Point - Triangle

/// @brief Compute the barycentric coordinates of the closest point on the triangle.
/// @param p Point
/// @param t0 Triangle's first vertex
/// @param t1 Triangle's second vertex
/// @param t2 Triangle's third vertex
/// @return Barycentric coordinates of the closest point
Eigen::Vector2d point_triangle_closest_point(
    const Eigen::Ref<const Eigen::Vector3d>& p,
    const Eigen::Ref<const Eigen::Vector3d>& t0,
    const Eigen::Ref<const Eigen::Vector3d>& t1,
    const Eigen::Ref<const Eigen::Vector3d>& t2);

/// @brief Compute the Jacobian of the closest point on the triangle.
/// @param p Point
/// @param t0 Triangle's first vertex
/// @param t1 Triangle's second vertex
/// @param t2 Triangle's third vertex
/// @return Jacobian of the closest point
Eigen::Matrix<double, 2, 12> point_triangle_closest_point_jacobian(
    const Eigen::Ref<const Eigen::Vector3d>& p,
    const Eigen::Ref<const Eigen::Vector3d>& t0,
    const Eigen::Ref<const Eigen::Vector3d>& t1,
    const Eigen::Ref<const Eigen::Vector3d>& t2);
// ============================================================================

namespace autogen {
    // J is (2×12) flattened in column-major order
    void edge_edge_closest_point_jacobian(
        double ea0_x,
        double ea0_y,
        double ea0_z,
        double ea1_x,
        double ea1_y,
        double ea1_z,
        double eb0_x,
        double eb0_y,
        double eb0_z,
        double eb1_x,
        double eb1_y,
        double eb1_z,
        double J[24]);

    // J is (2×12) flattened in column-major order
    void point_triangle_closest_point_jacobian(
        double p_x,
        double p_y,
        double p_z,
        double t0_x,
        double t0_y,
        double t0_z,
        double t1_x,
        double t1_y,
        double t1_z,
        double t2_x,
        double t2_y,
        double t2_z,
        double J[24]);
} // namespace autogen

} // namespace ipc
