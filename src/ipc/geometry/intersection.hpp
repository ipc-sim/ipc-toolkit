#pragma once

#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

/// @brief Check if an edge intersects a triangle.
/// @param e0 Edge start point.
/// @param e1 Edge end point.
/// @param t0 Triangle vertex 0.
/// @param t1 Triangle vertex 1.
/// @param t2 Triangle vertex 2.
/// @return True if the edge intersects the triangle.
bool is_edge_intersecting_triangle(
    Eigen::ConstRef<Eigen::Vector3d> e0,
    Eigen::ConstRef<Eigen::Vector3d> e1,
    Eigen::ConstRef<Eigen::Vector3d> t0,
    Eigen::ConstRef<Eigen::Vector3d> t1,
    Eigen::ConstRef<Eigen::Vector3d> t2);

} // namespace ipc
