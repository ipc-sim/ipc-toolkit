#pragma once

#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

// ============================================================================
// Point - Edge

double point_edge_closest_point(
    const Eigen::Ref<const VectorMax3d>& p,
    const Eigen::Ref<const VectorMax3d>& e0,
    const Eigen::Ref<const VectorMax3d>& e1);

VectorMax9d point_edge_closest_point_jacobian(
    const Eigen::Ref<const VectorMax3d>& p,
    const Eigen::Ref<const VectorMax3d>& e0,
    const Eigen::Ref<const VectorMax3d>& e1);

// ============================================================================
// Edge - Edge

/// Compute the barycentric coordinates of the closest points
Eigen::Vector2d edge_edge_closest_point(
    const Eigen::Ref<const Eigen::Vector3d>& ea0,
    const Eigen::Ref<const Eigen::Vector3d>& ea1,
    const Eigen::Ref<const Eigen::Vector3d>& eb0,
    const Eigen::Ref<const Eigen::Vector3d>& eb1);

Eigen::Matrix<double, 2, 12> edge_edge_closest_point_jacobian(
    const Eigen::Ref<const Eigen::Vector3d>& ea0,
    const Eigen::Ref<const Eigen::Vector3d>& ea1,
    const Eigen::Ref<const Eigen::Vector3d>& eb0,
    const Eigen::Ref<const Eigen::Vector3d>& eb1);

// ============================================================================
// Point - Triangle

/// Compute the barycentric coordinates of the closest point on the triangle.
Eigen::Vector2d point_triangle_closest_point(
    const Eigen::Ref<const Eigen::Vector3d>& p,
    const Eigen::Ref<const Eigen::Vector3d>& t0,
    const Eigen::Ref<const Eigen::Vector3d>& t1,
    const Eigen::Ref<const Eigen::Vector3d>& t2);

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
