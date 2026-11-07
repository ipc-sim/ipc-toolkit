#pragma once

#include <ipc/distance/distance_type.hpp>
#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/line_line.hpp>
#include <ipc/distance/point_line.hpp>
#include <ipc/math/math.hpp>
#include <ipc/utils/autodiff_types.hpp>
#include <ipc/utils/eigen_ext.hpp>
#include <ipc/utils/logger.hpp>

#include <algorithm>

namespace ipc {

/// @brief Edge-edge squared distance for an EdgeEdgeDistanceType that may
/// come from edge_edge_distance_type_exact.
/// @note The exact classifier has no parallel threshold, so it can return
/// EA_EB for parallel edges, where the line-line distance divides by zero.
/// Those cases fall back to the minimum point-line distance.
inline double edge_edge_distance_parallel_safe(
    Eigen::ConstRef<Eigen::Vector3d> ea0,
    Eigen::ConstRef<Eigen::Vector3d> ea1,
    Eigen::ConstRef<Eigen::Vector3d> eb0,
    Eigen::ConstRef<Eigen::Vector3d> eb1,
    const EdgeEdgeDistanceType dtype)
{
    if (dtype == EdgeEdgeDistanceType::EA_EB) {
        const Eigen::Vector3d normal = (ea1 - ea0).cross(eb1 - eb0);
        if (normal.squaredNorm() > 1e-20) {
            return line_line_distance(ea0, ea1, eb0, eb1);
        }
        return std::min(
            { point_line_distance(eb0, ea0, ea1),
              point_line_distance(eb1, ea0, ea1),
              point_line_distance(ea0, eb0, eb1),
              point_line_distance(ea1, eb0, eb1) });
    }
    return edge_edge_distance(ea0, ea1, eb0, eb1, dtype);
}

template <typename T>
T line_line_sqr_distance(
    Eigen::ConstRef<Eigen::Vector3<T>> ea0,
    Eigen::ConstRef<Eigen::Vector3<T>> ea1,
    Eigen::ConstRef<Eigen::Vector3<T>> eb0,
    Eigen::ConstRef<Eigen::Vector3<T>> eb1)
{
    const Eigen::Vector3<T> normal = (ea1 - ea0).cross(eb1 - eb0);
    const T line_to_line = (eb0 - ea0).dot(normal);
    return line_to_line * line_to_line / normal.squaredNorm();
}

template <typename T>
Eigen::Vector<T, 2> line_line_closest_point_pairs_uv(
    Eigen::ConstRef<Eigen::Vector3<T>> ea0,
    Eigen::ConstRef<Eigen::Vector3<T>> ea1,
    Eigen::ConstRef<Eigen::Vector3<T>> eb0,
    Eigen::ConstRef<Eigen::Vector3<T>> eb1)
{
    const Eigen::Vector3<T> u = ea1 - ea0;
    const Eigen::Vector3<T> v = eb1 - eb0;
    const Eigen::Vector3<T> w = ea0 - eb0;

    // We write this solve with cross products instead of the textbook dot
    // products, whose denominator is (u·u)(v·v) − (u·v)². Lagrange's identity
    // makes the two forms equal, but for nearly parallel edges (u·u)(v·v) and
    // (u·v)² agree to almost every digit, so their difference can round to
    // zero or below and the numerators lose all their digits too. That put
    // closest points edge lengths off the edge. The cross product u × v stays
    // accurate down to the builder's parallel threshold, which tests this same
    // ‖u × v‖², so every pair the builder keeps has a positive denominator
    // here.
    const Eigen::Vector3<T> n = u.cross(v);
    const T fac = n.squaredNorm();
    const T sN = n.dot(v.cross(w));
    const T tN = n.dot(u.cross(w));
    assert(fac > 0);

    return Eigen::Vector<T, 2>(sN, tN) / fac;
}

/// @brief Local coordinate on edge (ea0, ea1), clamped to [0, 1], of its point
/// closest to the line through (eb0, eb1). ESP edge-edge terms exist only for
/// EA_EB pairs, where this is the closest point between the edges.
template <typename T>
T closest_point_uv(
    Eigen::ConstRef<Eigen::Vector3<T>> ea0,
    Eigen::ConstRef<Eigen::Vector3<T>> ea1,
    Eigen::ConstRef<Eigen::Vector3<T>> eb0,
    Eigen::ConstRef<Eigen::Vector3<T>> eb1)
{
    T uv = line_line_closest_point_pairs_uv<T>(ea0, ea1, eb0, eb1)(0);

    if (uv < 0.) {
        uv = 0.;
    } else if (uv > 1.) {
        uv = 1.;
    }

    return uv;
}

} // namespace ipc
