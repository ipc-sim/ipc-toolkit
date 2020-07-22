#pragma once

#include <distance/distance_type.hpp>
#include <distance/line_line.hpp>
#include <distance/point_edge.hpp>
#include <distance/point_point.hpp>

namespace ipc {

/// @brief Compute the distance between a two lines segments in 3D.
/// @note The distance is actually squared distance.
/// @param ea0,ea1 The points of the first edge.
/// @param eb0,eb1 The points of the second edge.
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

    assert(false);
    return point_point_distance(ea0, eb0);
}

/// @brief Compute a mollifier for the edge-edge distance.
///
/// This helps smooth the non-smoothness at close to parallel edges.
///
/// @param ea0,ea1 The points of the first edge.
/// @param eb0,eb1 The points of the second edge.
template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1>
inline auto edge_edge_mollifier(
    const Eigen::MatrixBase<DerivedEA0>& ea0,
    const Eigen::MatrixBase<DerivedEA1>& ea1,
    const Eigen::MatrixBase<DerivedEB0>& eb0,
    const Eigen::MatrixBase<DerivedEB1>& eb1,
    const double eps_x)
{
    auto ee_cross_norm_sqr = (ea1 - ea0).cross(eb1 - eb0).squaredNorm();
    if (ee_cross_norm_sqr < eps_x) {
        auto x_div_eps_x = ee_cross_norm_sqr / eps_x;
        return (-x_div_eps_x + 2.0) * x_div_eps_x;
    } else {
        return decltype(ee_cross_norm_sqr)(1.0);
    }
}

/// @brief Compute the threshold of the mollifier edge-edge distance.
///
/// This values is computed based on the edges at rest length.
///
/// @param ea0,ea1 The points of the first edge at rest.
/// @param eb0,eb1 The points of the second edge at rest.
template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1>
double edge_edge_mollifier_threshold(
    const Eigen::MatrixBase<DerivedEA0>& ea0_rest,
    const Eigen::MatrixBase<DerivedEA1>& ea1_rest,
    const Eigen::MatrixBase<DerivedEB0>& eb0_rest,
    const Eigen::MatrixBase<DerivedEB1>& eb1_rest)
{
    return 1.0e-3 * (ea0_rest - ea1_rest).squaredNorm()
        * (eb0_rest - eb1_rest).squaredNorm();
}

} // namespace ipc
