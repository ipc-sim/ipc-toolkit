#pragma once

#include <distance/distance_type.hpp>
#include <distance/point_line.hpp>
#include <distance/point_point.hpp>

namespace ipc {

// http://geomalgorithms.com/a02-_lines.html
template <typename DerivedP, typename DerivedE0, typename DerivedE1>
auto point_edge_distance(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedE0>& e0,
    const Eigen::MatrixBase<DerivedE1>& e1)
{
    switch (point_edge_distance_type(p, e0, e1)) {
    case PointEdgeDistanceType::P_E0:
        return point_point_distance(p, e0);
    case PointEdgeDistanceType::P_E1:
        return point_point_distance(p, e1);
    case PointEdgeDistanceType::P_E:
        return point_line_distance(p, e0, e1);
    }
}

} // namespace ipc
