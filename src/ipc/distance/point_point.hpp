#pragma once

#include <ipc/config.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

/// @brief Compute the distance between two points.
/// @note The distance is actually squared distance.
/// @param p0 The first point.
/// @param p1 The second point.
/// @return The distance between p0 and p1.
IPC_TOOLKIT_HOST_DEVICE double point_point_distance(
    Eigen::ConstRef<VectorMax3d> p0, Eigen::ConstRef<VectorMax3d> p1);

/// @brief Compute the gradient of the distance between two points.
/// @note The distance is actually squared distance.
/// @param p0 The first point.
/// @param p1 The second point.
/// @return The computed gradient.
IPC_TOOLKIT_HOST_DEVICE VectorMax6d point_point_distance_gradient(
    Eigen::ConstRef<VectorMax3d> p0, Eigen::ConstRef<VectorMax3d> p1);

/// @brief Compute the hessian of the distance between two points.
/// @note The distance is actually squared distance.
/// @param p0 The first point.
/// @param p1 The second point.
/// @return The computed hessian.
IPC_TOOLKIT_HOST_DEVICE MatrixMax6d point_point_distance_hessian(
    Eigen::ConstRef<VectorMax3d> p0, Eigen::ConstRef<VectorMax3d> p1);

} // namespace ipc
