#pragma once

#include <distance/distance_type.hpp>
#include <distance/line_line.hpp>
#include <distance/point_edge.hpp>
#include <distance/point_point.hpp>

namespace ipc {

template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1>
auto edge_edge_distance(
    const Eigen::MatrixBase<DerivedEA0>& ea0,
    const Eigen::MatrixBase<DerivedEA1>& ea1,
    const Eigen::MatrixBase<DerivedEB0>& eb0,
    const Eigen::MatrixBase<DerivedEB1>& eb1)
{
    switch (edge_edge_distance_type(ea0, ea1, eb0, eb1)) {
    case EdgeEdgeDistanceType::EA0_EB0:
        return point_point_distance(ea0, eb0);

    case EdgeEdgeDistanceType::EA0_EB1:
        return point_point_distance(ea0, eb1);

    case EdgeEdgeDistanceType::EA1_EB0:
        return point_point_distance(ea1, eb0);

    case EdgeEdgeDistanceType::EA1_EB1:
        return point_point_distance(ea1, eb1);

    case EdgeEdgeDistanceType::EA_EB0:
        return point_edge_distance(eb0, ea0, ea1);

    case EdgeEdgeDistanceType::EA_EB1:
        return point_edge_distance(eb1, ea0, ea1);

    case EdgeEdgeDistanceType::EA0_EB:
        return point_edge_distance(ea0, eb0, eb1);

    case EdgeEdgeDistanceType::EA1_EB:
        return point_edge_distance(ea1, eb0, eb1);

    case EdgeEdgeDistanceType::EA_EB:
        return line_line_distance(ea0, ea1, eb0, eb1);
    }
}

} // namespace ipc
