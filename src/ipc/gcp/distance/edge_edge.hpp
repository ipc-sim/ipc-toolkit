#pragma once

#include "point_edge.hpp"

namespace ipc {

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

template <typename scalar>
scalar edge_edge_sqr_distance(
    Eigen::ConstRef<Eigen::Vector3<scalar>> ea0,
    Eigen::ConstRef<Eigen::Vector3<scalar>> ea1,
    Eigen::ConstRef<Eigen::Vector3<scalar>> eb0,
    Eigen::ConstRef<Eigen::Vector3<scalar>> eb1,
    EdgeEdgeDistanceType dtype)
{
    if constexpr (std::is_same<double, scalar>::value) {
        if (dtype == EdgeEdgeDistanceType::AUTO) {
            dtype = edge_edge_distance_type(ea0, ea1, eb0, eb1);
        }
    }

    switch (dtype) {
    case EdgeEdgeDistanceType::EA0_EB0:
        return PointEdgeDistance<scalar, 3>::point_point_sqr_distance(ea0, eb0);

    case EdgeEdgeDistanceType::EA0_EB1:
        return PointEdgeDistance<scalar, 3>::point_point_sqr_distance(ea0, eb1);

    case EdgeEdgeDistanceType::EA1_EB0:
        return PointEdgeDistance<scalar, 3>::point_point_sqr_distance(ea1, eb0);

    case EdgeEdgeDistanceType::EA1_EB1:
        return PointEdgeDistance<scalar, 3>::point_point_sqr_distance(ea1, eb1);

    case EdgeEdgeDistanceType::EA_EB0:
        return PointEdgeDistance<scalar, 3>::point_line_sqr_distance(
            eb0, ea0, ea1);

    case EdgeEdgeDistanceType::EA_EB1:
        return PointEdgeDistance<scalar, 3>::point_line_sqr_distance(
            eb1, ea0, ea1);

    case EdgeEdgeDistanceType::EA0_EB:
        return PointEdgeDistance<scalar, 3>::point_line_sqr_distance(
            ea0, eb0, eb1);

    case EdgeEdgeDistanceType::EA1_EB:
        return PointEdgeDistance<scalar, 3>::point_line_sqr_distance(
            ea1, eb0, eb1);

    case EdgeEdgeDistanceType::EA_EB:
        return line_line_sqr_distance<scalar>(ea0, ea1, eb0, eb1);

    default:
        throw std::invalid_argument(
            "Invalid distance type for edge-edge distance!");
    }
}

template <typename T>
Eigen::Vector3<T> line_line_closest_point_direction(
    Eigen::ConstRef<Eigen::Vector3<T>> ea0,
    Eigen::ConstRef<Eigen::Vector3<T>> ea1,
    Eigen::ConstRef<Eigen::Vector3<T>> eb0,
    Eigen::ConstRef<Eigen::Vector3<T>> eb1);

std::tuple<Eigen::Vector3d, Eigen::Matrix<double, 3, 12>>
line_line_closest_point_direction_gradient(
    Eigen::ConstRef<Eigen::Vector3d> ea0,
    Eigen::ConstRef<Eigen::Vector3d> ea1,
    Eigen::ConstRef<Eigen::Vector3d> eb0,
    Eigen::ConstRef<Eigen::Vector3d> eb1);

std::tuple<
    Eigen::Vector3d,
    Eigen::Matrix<double, 3, 12>,
    std::array<Matrix12d, 3>>
line_line_closest_point_direction_hessian(
    Eigen::ConstRef<Eigen::Vector3d> ea0,
    Eigen::ConstRef<Eigen::Vector3d> ea1,
    Eigen::ConstRef<Eigen::Vector3d> eb0,
    Eigen::ConstRef<Eigen::Vector3d> eb1);

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

    const T a = u.squaredNorm();
    const T b = u.dot(v);
    const T c = v.squaredNorm();
    const T d = u.dot(w);
    const T e = v.dot(w);

    const T sN = (b * e - c * d);
    const T tN = (a * e - b * d);
    const T fac = a * c - pow(b, 2);
    assert(fac > 0);

    return Eigen::Vector<T, 2>(sN, tN) / fac;
}

template <typename T>
Eigen::Matrix<T, 3, 2> line_line_closest_point_pairs(
    Eigen::ConstRef<Eigen::Vector3<T>> ea0,
    Eigen::ConstRef<Eigen::Vector3<T>> ea1,
    Eigen::ConstRef<Eigen::Vector3<T>> eb0,
    Eigen::ConstRef<Eigen::Vector3<T>> eb1)
{
    const Eigen::Vector<T, 2> uvs =
        line_line_closest_point_pairs_uv<T>(ea0, ea1, eb0, eb1);

    Eigen::Matrix<T, 3, 2> out;
    out.col(0) = ea0 + uvs(0) * (ea1 - ea0);
    out.col(1) = eb0 + uvs(1) * (eb1 - eb0);

    return out;
}

std::tuple<Vector6d, Eigen::Matrix<double, 6, 12>>
line_line_closest_point_pairs_gradient(
    Eigen::ConstRef<Eigen::Vector3d> ea0,
    Eigen::ConstRef<Eigen::Vector3d> ea1,
    Eigen::ConstRef<Eigen::Vector3d> eb0,
    Eigen::ConstRef<Eigen::Vector3d> eb1);

std::tuple<Vector6d, Eigen::Matrix<double, 6, 12>, std::array<Matrix12d, 6>>
line_line_closest_point_pairs_hessian(
    Eigen::ConstRef<Eigen::Vector3d> ea0,
    Eigen::ConstRef<Eigen::Vector3d> ea1,
    Eigen::ConstRef<Eigen::Vector3d> eb0,
    Eigen::ConstRef<Eigen::Vector3d> eb1);

/// @brief Computes the direction of the closest point pair
/// @param ea0 Vertex 0 of edge 0
/// @param ea1 Vertex 1 of edge 0
/// @param eb0 Vertex 0 of edge 1
/// @param eb1 Vertex 1 of edge 1
/// @param dtype Edge-edge distance type
/// @return Difference of the pair of closest point, pointing from edge 0 to edge 1
template <typename T>
Eigen::Vector3<T> edge_edge_closest_point_direction(
    Eigen::ConstRef<Eigen::Vector3<T>> ea0,
    Eigen::ConstRef<Eigen::Vector3<T>> ea1,
    Eigen::ConstRef<Eigen::Vector3<T>> eb0,
    Eigen::ConstRef<Eigen::Vector3<T>> eb1,
    EdgeEdgeDistanceType dtype);

/// @brief Computes the position of two closest points on two edges
/// @param ea0 Vertex 0 of edge 0
/// @param ea1 Vertex 1 of edge 0
/// @param eb0 Vertex 0 of edge 1
/// @param eb1 Vertex 1 of edge 1
/// @param dtype Edge-edge distance type
template <typename T>
Eigen::Matrix<T, 3, 2> edge_edge_closest_point_pairs(
    Eigen::ConstRef<Eigen::Vector3<T>> ea0,
    Eigen::ConstRef<Eigen::Vector3<T>> ea1,
    Eigen::ConstRef<Eigen::Vector3<T>> eb0,
    Eigen::ConstRef<Eigen::Vector3<T>> eb1,
    EdgeEdgeDistanceType dtype);

// Compute the closest point local coordinate on edge (e0, e1) with respect to
// edge (e2, e3) This function is written in a consistent way as the edge-edge
// distance type classification
template <typename T>
T closest_point_uv(
    Eigen::ConstRef<Eigen::Vector3<T>> e0,
    Eigen::ConstRef<Eigen::Vector3<T>> e1,
    Eigen::ConstRef<Eigen::Vector3<T>> e2,
    Eigen::ConstRef<Eigen::Vector3<T>> e3,
    EdgeEdgeDistanceType dtype)
{
    Eigen::Vector<T, 3> u = e1 - e0;
    Eigen::Vector<T, 3> v = e3 - e2;

    T uv(0.);
    if (dtype == EdgeEdgeDistanceType::EA_EB) {
        Eigen::Vector2<T> uvs =
            line_line_closest_point_pairs_uv<T>(e0, e1, e2, e3);

        uv = uvs(0);
    } else if (dtype == EdgeEdgeDistanceType::EA_EB0) {
        const T a = u.squaredNorm();
        const T d = u.dot(e0 - e2);
        uv = (-d) / a;
    } else if (dtype == EdgeEdgeDistanceType::EA_EB1) {
        const T a = u.squaredNorm();
        const T b = u.dot(v);
        const T d = u.dot(e0 - e2);
        uv = (-d + b) / a;
    } else
        log_and_throw_error(
            "edge-edge dtype {} cannot handle!", static_cast<int>(dtype));

    if (uv < 0.) {
        uv = 0.;
    } else if (uv > 1.) {
        uv = 1.;
    }

    return uv;
}
} // namespace ipc
