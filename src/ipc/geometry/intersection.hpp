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

/// @brief Edge–triangle intersection test that also returns the hit location.
///
/// Uses the same robust orient3d plane-side gate as
/// `is_edge_intersecting_triangle`, then solves for the barycentric (u, v) on
/// the triangle and the parameter t along the edge. Boundary hits count as
/// intersections (the comparisons are inclusive), matching
/// `is_edge_intersecting_triangle`.
///
/// @note The out-parameters are only meaningful when this returns true, and
/// are set to NaN otherwise. The one exception is a degenerate configuration
/// (edge coplanar with the triangle, or a degenerate triangle), where this
/// conservatively returns true but the coordinates are not uniquely defined
/// and are left as NaN. Callers that need the hit point must check for NaN.
///
/// @param[in] e0,e1 Edge endpoints.
/// @param[in] t0,t1,t2 Triangle vertices.
/// @param[out] u,v Triangle barycentric coordinates (along t1−t0, t2−t0).
/// @param[out] t Edge parameter in [0, 1].
/// @return True if the edge intersects the triangle (including its boundary).
bool edge_triangle_intersection(
    Eigen::ConstRef<Eigen::Vector3d> e0,
    Eigen::ConstRef<Eigen::Vector3d> e1,
    Eigen::ConstRef<Eigen::Vector3d> t0,
    Eigen::ConstRef<Eigen::Vector3d> t1,
    Eigen::ConstRef<Eigen::Vector3d> t2,
    double& u,
    double& v,
    double& t);

} // namespace ipc
