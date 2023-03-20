#pragma once

#include <ipc/distance/distance_type.hpp>

namespace ipc {

/// @brief Compute the distance between a two lines segments in 3D.
/// @note The distance is actually squared distance.
/// @param ea0 The first vertex of the first edge.
/// @param ea1 The second vertex of the first edge.
/// @param eb0 The first vertex of the second edge.
/// @param eb1 The second vertex of the second edge.
/// @param dtype The point edge distance type to compute.
/// @return The distance between the two edges.
double edge_edge_distance(
    const Eigen::Ref<const Eigen::Vector3d>& ea0,
    const Eigen::Ref<const Eigen::Vector3d>& ea1,
    const Eigen::Ref<const Eigen::Vector3d>& eb0,
    const Eigen::Ref<const Eigen::Vector3d>& eb1,
    EdgeEdgeDistanceType dtype = EdgeEdgeDistanceType::AUTO);

/// @brief Compute the gradient of the distance between a two lines segments.
/// @note The distance is actually squared distance.
/// @param ea0 The first vertex of the first edge.
/// @param ea1 The second vertex of the first edge.
/// @param eb0 The first vertex of the second edge.
/// @param eb1 The second vertex of the second edge.
/// @param dtype The point edge distance type to compute.
/// @return The gradient of the distance wrt ea0, ea1, eb0, and eb1.
Vector12d edge_edge_distance_gradient(
    const Eigen::Ref<const Eigen::Vector3d>& ea0,
    const Eigen::Ref<const Eigen::Vector3d>& ea1,
    const Eigen::Ref<const Eigen::Vector3d>& eb0,
    const Eigen::Ref<const Eigen::Vector3d>& eb1,
    EdgeEdgeDistanceType dtype = EdgeEdgeDistanceType::AUTO);

/// @brief Compute the hessian of the distance between a two lines segments.
/// @note The distance is actually squared distance.
/// @param ea0 The first vertex of the first edge.
/// @param ea1 The second vertex of the first edge.
/// @param eb0 The first vertex of the second edge.
/// @param eb1 The second vertex of the second edge.
/// @param dtype The point edge distance type to compute.
/// @return The hessian of the distance wrt ea0, ea1, eb0, and eb1.
Matrix12d edge_edge_distance_hessian(
    const Eigen::Ref<const Eigen::Vector3d>& ea0,
    const Eigen::Ref<const Eigen::Vector3d>& ea1,
    const Eigen::Ref<const Eigen::Vector3d>& eb0,
    const Eigen::Ref<const Eigen::Vector3d>& eb1,
    EdgeEdgeDistanceType dtype = EdgeEdgeDistanceType::AUTO);

} // namespace ipc
