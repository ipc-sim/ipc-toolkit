#pragma once

#include <distance/distance_type.hpp>
#include <distance/point_line.hpp>
#include <distance/point_plane.hpp>
#include <distance/point_point.hpp>

namespace ipc {

template <
    typename DerivedP,
    typename DerivedT0,
    typename DerivedT1,
    typename DerivedT2>
auto point_triangle_distance(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedT0>& t0,
    const Eigen::MatrixBase<DerivedT1>& t1,
    const Eigen::MatrixBase<DerivedT2>& t2)
{
    switch (point_triangle_distance_type(p, t0, t1, t2)) {
    case PointTriangleDistanceType::P_T0:
        return point_point_distance(p, t0);

    case PointTriangleDistanceType::P_T1:
        return point_point_distance(p, t1);

    case PointTriangleDistanceType::P_T2:
        return point_point_distance(p, t2);

    case PointTriangleDistanceType::P_E0:
        return point_line_distance(p, t0, t1);

    case PointTriangleDistanceType::P_E1:
        return point_line_distance(p, t1, t2);

    case PointTriangleDistanceType::P_E2:
        return point_line_distance(p, t2, t0);

    case PointTriangleDistanceType::P_T:
        return point_plane_distance(p, t0, t1, t2);
    }
}

} // namespace ipc
