#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>

namespace ipc {

// NOTE: squared distance

// Compute the distance between a two lines and a plane (defined by a triangle)
template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1>
auto line_line_distance(
    const Eigen::MatrixBase<DerivedEA0>& ea0,
    const Eigen::MatrixBase<DerivedEA1>& ea1,
    const Eigen::MatrixBase<DerivedEB0>& eb0,
    const Eigen::MatrixBase<DerivedEB1>& eb1)
{
    const auto normal = (ea1 - ea0).cross(eb1 - eb0);
    if (normal.squaredNorm() == 0) {
        throw "parallel lines";
    }
#ifdef USE_DISTANCE_SQUARED
    const auto line_to_line = (eb0 - ea0).dot(normal);
    return line_to_line * line_to_line / normal.squaredNorm();
#else
    normal.normalize();
    return (eb0 - ea0).dot(normal);
#endif
}

} // namespace ipc
