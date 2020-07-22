#pragma once

#include <distance/distance_type.hpp>
#include <distance/point_line.hpp>
#include <distance/point_point.hpp>

namespace ipc {

/// @brief Compute the distance between a point and edge in 2D or 3D.
/// @note The distance is actually squared distance.
/// @param p The point.
/// @param e0,e1 The points of the edge.
template <typename DerivedP, typename DerivedE0, typename DerivedE1>
auto point_edge_distance(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedE0>& e0,
    const Eigen::MatrixBase<DerivedE1>& e1)
{
    // switch (point_edge_distance_type(p, e0, e1)) {
    // case PointEdgeDistanceType::P_E0:
    //     return point_point_distance(p, e0);
    // case PointEdgeDistanceType::P_E1:
    //     return point_point_distance(p, e1);
    // case PointEdgeDistanceType::P_E:
    //     return point_line_distance(p, e0, e1);
    // }
    typedef typename DerivedP::Scalar T;

    // Project the point onto the line
    auto e = e1 - e0;
    auto e_length_sqr = e.squaredNorm();
    auto alpha = e_length_sqr != 0 ? ((p - e0).dot(e) / e_length_sqr) : T(0.5);

    if (alpha <= 0) {
        return point_point_distance(p, e0);
    } else if (alpha >= 1) {
        return point_point_distance(p, e1);
    }
    return point_point_distance(p, e * alpha + e0);
}

} // namespace ipc
