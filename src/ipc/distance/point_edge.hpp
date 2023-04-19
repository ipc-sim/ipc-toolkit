#pragma once

#include <ipc/distance/distance_type.hpp>

namespace ipc {

/// @brief Compute the distance between a point and edge in 2D or 3D.
/// @note The distance is actually squared distance.
/// @param p The point.
/// @param e0 The first vertex of the edge.
/// @param e1 The second vertex of the edge.
/// @param dtype The point edge distance type to compute.
/// @return The distance between the point and edge.
double point_edge_distance(
    const Eigen::Ref<const VectorMax3d>& p,
    const Eigen::Ref<const VectorMax3d>& e0,
    const Eigen::Ref<const VectorMax3d>& e1,
    PointEdgeDistanceType dtype = PointEdgeDistanceType::AUTO);

/// @brief Compute the gradient of the distance between a point and edge.
/// @note The distance is actually squared distance.
/// @param p The point.
/// @param e0 The first vertex of the edge.
/// @param e1 The second vertex of the edge.
/// @param dtype The point edge distance type to compute.
/// @return grad The gradient of the distance wrt p, e0, and e1.
VectorMax9d point_edge_distance_gradient(
    const Eigen::Ref<const VectorMax3d>& p,
    const Eigen::Ref<const VectorMax3d>& e0,
    const Eigen::Ref<const VectorMax3d>& e1,
    PointEdgeDistanceType dtype = PointEdgeDistanceType::AUTO);

/// @brief Compute the hessian of the distance between a point and edge.
/// @note The distance is actually squared distance.
/// @param p The point.
/// @param e0 The first vertex of the edge.
/// @param e1 The second vertex of the edge.
/// @param dtype The point edge distance type to compute.
/// @return hess The hessian of the distance wrt p, e0, and e1.
MatrixMax9d point_edge_distance_hessian(
    const Eigen::Ref<const VectorMax3d>& p,
    const Eigen::Ref<const VectorMax3d>& e0,
    const Eigen::Ref<const VectorMax3d>& e1,
    PointEdgeDistanceType dtype = PointEdgeDistanceType::AUTO);

} // namespace ipc
