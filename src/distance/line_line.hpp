#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>

namespace ipc {

/// @brief Compute the distance between a two infinite lines in 3D.
/// @note The distance is actually squared distance.
/// @warning If the lines are parallel this function returns a distance of zero.
/// @param ea0,ea1 The points of the edge defining the first line.
/// @param eb0,eb1 The points of the edge defining the second line.
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
    const auto line_to_line = (eb0 - ea0).dot(normal);
    return line_to_line * line_to_line / normal.squaredNorm();
}

} // namespace ipc
