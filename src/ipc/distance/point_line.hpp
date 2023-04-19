#pragma once

#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

/// @brief Compute the distance between a point and line in 2D or 3D.
/// @note The distance is actually squared distance.
/// @param p The point.
/// @param e0 The first vertex of the edge defining the line.
/// @param e1 The second vertex of the edge defining the line.
/// @return The distance between the point and line.
double point_line_distance(
    const Eigen::Ref<const VectorMax3d>& p,
    const Eigen::Ref<const VectorMax3d>& e0,
    const Eigen::Ref<const VectorMax3d>& e1);

/// @brief Compute the gradient of the distance between a point and line.
/// @note The distance is actually squared distance.
/// @param p The point.
/// @param e0 The first vertex of the edge defining the line.
/// @param e1 The second vertex of the edge defining the line.
/// @return The gradient of the distance wrt p, e0, and e1.
VectorMax9d point_line_distance_gradient(
    const Eigen::Ref<const VectorMax3d>& p,
    const Eigen::Ref<const VectorMax3d>& e0,
    const Eigen::Ref<const VectorMax3d>& e1);

/// @brief Compute the hessian of the distance between a point and line.
/// @note The distance is actually squared distance.
/// @param p The point.
/// @param e0 The first vertex of the edge defining the line.
/// @param e1 The second vertex of the edge defining the line.
/// @return The hessian of the distance wrt p, e0, and e1.
MatrixMax9d point_line_distance_hessian(
    const Eigen::Ref<const VectorMax3d>& p,
    const Eigen::Ref<const VectorMax3d>& e0,
    const Eigen::Ref<const VectorMax3d>& e1);

// Symbolically generated derivatives;
namespace autogen {
    void point_line_distance_gradient_2D(
        double v01,
        double v02,
        double v11,
        double v12,
        double v21,
        double v22,
        double g[6]);

    void point_line_distance_gradient_3D(
        double v01,
        double v02,
        double v03,
        double v11,
        double v12,
        double v13,
        double v21,
        double v22,
        double v23,
        double g[9]);

    void point_line_distance_hessian_2D(
        double v01,
        double v02,
        double v11,
        double v12,
        double v21,
        double v22,
        double H[36]);

    void point_line_distance_hessian_3D(
        double v01,
        double v02,
        double v03,
        double v11,
        double v12,
        double v13,
        double v21,
        double v22,
        double v23,
        double H[81]);
} // namespace autogen
} // namespace ipc
