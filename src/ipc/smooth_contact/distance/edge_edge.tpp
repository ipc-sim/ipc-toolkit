#pragma once

#include "edge_edge.hpp"

namespace ipc {

template <typename scalar>
scalar line_line_sqr_distance(
    Eigen::ConstRef<Vector3<scalar>> ea0,
    Eigen::ConstRef<Vector3<scalar>> ea1,
    Eigen::ConstRef<Vector3<scalar>> eb0,
    Eigen::ConstRef<Vector3<scalar>> eb1)
{
    const Vector3<scalar> normal = (ea1 - ea0).cross(eb1 - eb0);
    const scalar line_to_line = (eb0 - ea0).dot(normal);
    return line_to_line * line_to_line / normal.squaredNorm();
}

template <typename scalar>
scalar edge_edge_sqr_distance(
    Eigen::ConstRef<Vector3<scalar>> ea0,
    Eigen::ConstRef<Vector3<scalar>> ea1,
    Eigen::ConstRef<Vector3<scalar>> eb0,
    Eigen::ConstRef<Vector3<scalar>> eb1,
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

template <typename scalar>
Vector3<scalar> line_line_closest_point_direction(
    Eigen::ConstRef<Vector3<scalar>> ea0,
    Eigen::ConstRef<Vector3<scalar>> ea1,
    Eigen::ConstRef<Vector3<scalar>> eb0,
    Eigen::ConstRef<Vector3<scalar>> eb1)
{
    const Vector3<scalar> normal = (ea1 - ea0).cross(eb1 - eb0);
    return ((eb0 - ea0).dot(normal) / normal.squaredNorm()) * normal;
}

template <typename scalar>
Eigen::Matrix<scalar, 3, 2> line_line_closest_point_pairs(
    Eigen::ConstRef<Vector3<scalar>> ea0,
    Eigen::ConstRef<Vector3<scalar>> ea1,
    Eigen::ConstRef<Vector3<scalar>> eb0,
    Eigen::ConstRef<Vector3<scalar>> eb1)
{
    const Vector3<scalar> ta = ea1 - ea0;
    const Vector3<scalar> tb = eb1 - eb0;
    const scalar la = ta.squaredNorm();
    const scalar lb = tb.squaredNorm();
    const scalar lab = ta.dot(tb);
    const Vector3<scalar> d = eb0 - ea0;

    Eigen::Matrix<scalar, 3, 2> out;
    const scalar fac = la * lb - pow(lab, 2);
    out.col(0) = ea0 + (lb * ta.dot(d) - lab * tb.dot(d)) / fac * ta;
    out.col(1) = eb0 + (lab * ta.dot(d) - la * tb.dot(d)) / fac * tb;

    return out;
}

/// @brief Computes the direction of the closest point pair
/// @param ea0 Vertex 0 of edge 0
/// @param ea1 Vertex 1 of edge 0
/// @param eb0 Vertex 0 of edge 1
/// @param eb1 Vertex 1 of edge 1
/// @param dtype Edge-edge distance type
/// @return Difference of the pair of closest point, pointing from edge 0 to edge 1
template <typename scalar>
Vector3<scalar> edge_edge_closest_point_direction(
    Eigen::ConstRef<Vector3<scalar>> ea0,
    Eigen::ConstRef<Vector3<scalar>> ea1,
    Eigen::ConstRef<Vector3<scalar>> eb0,
    Eigen::ConstRef<Vector3<scalar>> eb1,
    EdgeEdgeDistanceType dtype)
{
    if constexpr (std::is_same<double, scalar>::value) {
        if (dtype == EdgeEdgeDistanceType::AUTO) {
            dtype = edge_edge_distance_type(ea0, ea1, eb0, eb1);
        }
    }

    switch (dtype) {
    case EdgeEdgeDistanceType::EA0_EB0:
        return (eb0 - ea0);

    case EdgeEdgeDistanceType::EA0_EB1:
        return (eb1 - ea0);

    case EdgeEdgeDistanceType::EA1_EB0:
        return (eb0 - ea1);

    case EdgeEdgeDistanceType::EA1_EB1:
        return (eb1 - ea1);

    case EdgeEdgeDistanceType::EA_EB0:
        return PointEdgeDistance<scalar, 3>::point_line_closest_point_direction(
            eb0, ea0, ea1);

    case EdgeEdgeDistanceType::EA_EB1:
        return PointEdgeDistance<scalar, 3>::point_line_closest_point_direction(
            eb1, ea0, ea1);

    case EdgeEdgeDistanceType::EA0_EB:
        return -PointEdgeDistance<
            scalar, 3>::point_line_closest_point_direction(ea0, eb0, eb1);

    case EdgeEdgeDistanceType::EA1_EB:
        return -PointEdgeDistance<
            scalar, 3>::point_line_closest_point_direction(ea1, eb0, eb1);

    case EdgeEdgeDistanceType::EA_EB:
        return line_line_closest_point_direction<scalar>(ea0, ea1, eb0, eb1);

    default:
        throw std::invalid_argument(
            "Invalid distance type for edge-edge distance!");
    }
}

template <typename scalar>
Eigen::Matrix<scalar, 3, 2> edge_edge_closest_point_pairs(
    Eigen::ConstRef<Vector3<scalar>> ea0,
    Eigen::ConstRef<Vector3<scalar>> ea1,
    Eigen::ConstRef<Vector3<scalar>> eb0,
    Eigen::ConstRef<Vector3<scalar>> eb1,
    EdgeEdgeDistanceType dtype)
{
    if constexpr (std::is_same<double, scalar>::value) {
        if (dtype == EdgeEdgeDistanceType::AUTO) {
            dtype = edge_edge_distance_type(ea0, ea1, eb0, eb1);
        }
    }

    Eigen::Matrix<scalar, 3, 2> out;
    switch (dtype) {
    case EdgeEdgeDistanceType::EA0_EB0:
        out << ea0, eb0;
        break;

    case EdgeEdgeDistanceType::EA0_EB1:
        out << ea0, eb1;
        break;

    case EdgeEdgeDistanceType::EA1_EB0:
        out << ea1, eb0;
        break;

    case EdgeEdgeDistanceType::EA1_EB1:
        out << ea1, eb1;
        break;

    case EdgeEdgeDistanceType::EA_EB0:
        out << eb0
                - PointEdgeDistance<scalar, 3>::
                    point_line_closest_point_direction(eb0, ea0, ea1),
            eb0;
        break;

    case EdgeEdgeDistanceType::EA_EB1:
        out << eb1
                - PointEdgeDistance<scalar, 3>::
                    point_line_closest_point_direction(eb1, ea0, ea1),
            eb1;
        break;

    case EdgeEdgeDistanceType::EA0_EB:
        out << ea0,
            ea0
            - PointEdgeDistance<scalar, 3>::point_line_closest_point_direction(
                ea0, eb0, eb1);
        break;

    case EdgeEdgeDistanceType::EA1_EB:
        out << ea1,
            ea1
            - PointEdgeDistance<scalar, 3>::point_line_closest_point_direction(
                ea1, eb0, eb1);
        break;

    case EdgeEdgeDistanceType::EA_EB:
        out = line_line_closest_point_pairs<scalar>(ea0, ea1, eb0, eb1);
        break;

    default:
        throw std::invalid_argument(
            "Invalid distance type for edge-edge distance!");
    }

    return out;
}

} // namespace ipc