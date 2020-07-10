#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>

namespace ipc {

// NOTE: squared distance

// Compute the distance between a point and a plane (defined by a triangle)
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
#ifdef USE_DISTANCE_SQUARED
    auto point_to_plane = (p - t0).dot(normal);
    return point_to_plane * point_to_plane / normal.squaredNorm();
#else
    return (p - t0).dot(normal.normalized());
#endif
}

} // namespace ipc
