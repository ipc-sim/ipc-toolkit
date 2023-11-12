#pragma once

#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

/// @brief Compute the distance between a two infinite lines in 3D.
/// @note The distance is actually squared distance.
/// @warning If the lines are parallel this function returns a distance of zero.
/// @param ea0 The first vertex of the edge defining the first line.
/// @param ea1 The second vertex of the edge defining the first line.
/// @param ea0 The first vertex of the edge defining the second line.
/// @param ea1 The second vertex of the edge defining the second line.
/// @return The distance between the two lines.
double line_line_distance(
    const Eigen::Ref<const Eigen::Vector3d>& ea0,
    const Eigen::Ref<const Eigen::Vector3d>& ea1,
    const Eigen::Ref<const Eigen::Vector3d>& eb0,
    const Eigen::Ref<const Eigen::Vector3d>& eb1);

/// @brief Compute the gradient of the distance between a two lines in 3D.
/// @note The distance is actually squared distance.
/// @warning If the lines are parallel this function returns a distance of zero.
/// @param ea0 The first vertex of the edge defining the first line.
/// @param ea1 The second vertex of the edge defining the first line.
/// @param ea0 The first vertex of the edge defining the second line.
/// @param ea1 The second vertex of the edge defining the second line.
/// @return The gradient of the distance wrt ea0, ea1, eb0, and eb1.
Vector12d line_line_distance_gradient(
    const Eigen::Ref<const Eigen::Vector3d>& ea0,
    const Eigen::Ref<const Eigen::Vector3d>& ea1,
    const Eigen::Ref<const Eigen::Vector3d>& eb0,
    const Eigen::Ref<const Eigen::Vector3d>& eb1);

/// @brief Compute the hessian of the distance between a two lines in 3D.
/// @note The distance is actually squared distance.
/// @warning If the lines are parallel this function returns a distance of zero.
/// @param ea0 The first vertex of the edge defining the first line.
/// @param ea1 The second vertex of the edge defining the first line.
/// @param ea0 The first vertex of the edge defining the second line.
/// @param ea1 The second vertex of the edge defining the second line.
/// @return The hessian of the distance wrt ea0, ea1, eb0, and eb1.
Matrix12d line_line_distance_hessian(
    const Eigen::Ref<const Eigen::Vector3d>& ea0,
    const Eigen::Ref<const Eigen::Vector3d>& ea1,
    const Eigen::Ref<const Eigen::Vector3d>& eb0,
    const Eigen::Ref<const Eigen::Vector3d>& eb1);

// Symbolically generated derivatives;
namespace autogen {
    void line_line_distance_gradient(
        double v01,
        double v02,
        double v03,
        double v11,
        double v12,
        double v13,
        double v21,
        double v22,
        double v23,
        double v31,
        double v32,
        double v33,
        double g[12]);

    void line_line_distance_hessian(
        double v01,
        double v02,
        double v03,
        double v11,
        double v12,
        double v13,
        double v21,
        double v22,
        double v23,
        double v31,
        double v32,
        double v33,
        double H[144]);
} // namespace autogen

} // namespace ipc
