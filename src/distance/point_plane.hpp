#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>

namespace ipc {

/// @brief Compute the distance between a point and a plane (defined by a
/// triangle).
/// @note The distance is actually squared distance.
/// @param p The point.
/// @param t0,t1,t2 The points of the triangle defining the plane.
template <
    typename DerivedP,
    typename DerivedT0,
    typename DerivedT1,
    typename DerivedT2>
auto point_plane_distance(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedT0>& t0,
    const Eigen::MatrixBase<DerivedT1>& t1,
    const Eigen::MatrixBase<DerivedT2>& t2)
{
    auto normal = (t1 - t0).cross(t2 - t0);
    auto point_to_plane = (p - t0).dot(normal);
    return point_to_plane * point_to_plane / normal.squaredNorm();
}

} // namespace ipc
