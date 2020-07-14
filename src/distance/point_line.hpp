#pragma once

#include <Eigen/Core>

#include <distance/point_point.hpp>

namespace ipc {

/// @brief Compute the distance between a point and a line (defined by an edge).
/// @note The distance is actually squared distance.
/// @param p The point.
/// @param e0,e1 The points of the edge defining the line.
template <typename DerivedP, typename DerivedE0, typename DerivedE1>
auto point_line_distance(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedE0>& e0,
    const Eigen::MatrixBase<DerivedE1>& e1)
{
    typedef typename DerivedP::Scalar T;

    // Project the point onto the line
    auto e = e1 - e0;
    auto e_length_sqr = e.squaredNorm();
    auto alpha = e_length_sqr != 0 ? ((p - e0).dot(e) / e_length_sqr) : T(0.5);

    return point_point_distance(p, e * alpha + e0);
}

} // namespace ipc
