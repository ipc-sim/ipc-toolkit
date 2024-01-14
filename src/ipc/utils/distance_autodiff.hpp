#pragma once

#include <ipc/distance/distance_type.hpp>
#include "math.hpp"

namespace ipc {
    template <typename scalar>
    scalar point_point_sqr_distance(
        const Eigen::Ref<const Vector3<scalar>>& a,
        const Eigen::Ref<const Vector3<scalar>>& b)
    {
        return (a - b).squaredNorm();
    }

    template <typename scalar>
    scalar line_line_sqr_distance(
        const Eigen::Ref<const Vector3<scalar>>& ea0,
        const Eigen::Ref<const Vector3<scalar>>& ea1,
        const Eigen::Ref<const Vector3<scalar>>& eb0,
        const Eigen::Ref<const Vector3<scalar>>& eb1)
    {
        const Vector3<scalar> normal = (ea1 - ea0).cross(eb1 - eb0);
        const scalar line_to_line = (eb0 - ea0).dot(normal);
        return line_to_line * line_to_line / normal.squaredNorm();
    }

    template <typename scalar>
    scalar point_line_sqr_distance(
        const Eigen::Ref<const Vector3<scalar>>& p,
        const Eigen::Ref<const Vector3<scalar>>& e0,
        const Eigen::Ref<const Vector3<scalar>>& e1)
    {
        return (e0 - p).cross(e1 - p).squaredNorm() / (e1 - e0).squaredNorm();
    }

    template <typename scalar>
    scalar edge_edge_sqr_distance(
        const Eigen::Ref<const Vector3<scalar>>& ea0,
        const Eigen::Ref<const Vector3<scalar>>& ea1,
        const Eigen::Ref<const Vector3<scalar>>& eb0,
        const Eigen::Ref<const Vector3<scalar>>& eb1,
        EdgeEdgeDistanceType dtype)
    {
        switch (dtype) {
        case EdgeEdgeDistanceType::EA0_EB0:
            return point_point_sqr_distance<scalar>(ea0, eb0);

        case EdgeEdgeDistanceType::EA0_EB1:
            return point_point_sqr_distance<scalar>(ea0, eb1);

        case EdgeEdgeDistanceType::EA1_EB0:
            return point_point_sqr_distance<scalar>(ea1, eb0);

        case EdgeEdgeDistanceType::EA1_EB1:
            return point_point_sqr_distance<scalar>(ea1, eb1);

        case EdgeEdgeDistanceType::EA_EB0:
            return point_line_sqr_distance<scalar>(eb0, ea0, ea1);

        case EdgeEdgeDistanceType::EA_EB1:
            return point_line_sqr_distance<scalar>(eb1, ea0, ea1);

        case EdgeEdgeDistanceType::EA0_EB:
            return point_line_sqr_distance<scalar>(ea0, eb0, eb1);

        case EdgeEdgeDistanceType::EA1_EB:
            return point_line_sqr_distance<scalar>(ea1, eb0, eb1);

        case EdgeEdgeDistanceType::EA_EB:
            return line_line_sqr_distance<scalar>(ea0, ea1, eb0, eb1);

        default:
            throw std::invalid_argument(
                "Invalid distance type for edge-edge distance!");
        }
    }

    template <typename scalar>
    Vector3<scalar> point_line_closest_point_direction(
        const Eigen::Ref<const Vector3<scalar>>& p,
        const Eigen::Ref<const Vector3<scalar>>& e0,
        const Eigen::Ref<const Vector3<scalar>>& e1)
    {
        const Vector3<scalar> tangent = (e1 - e0).normalized();
        return (p - e0) - (p - e0).dot(tangent) * tangent;
    }

    template <typename scalar>
    Vector3<scalar> line_line_closest_point_direction(
        const Eigen::Ref<const Vector3<scalar>>& ea0,
        const Eigen::Ref<const Vector3<scalar>>& ea1,
        const Eigen::Ref<const Vector3<scalar>>& eb0,
        const Eigen::Ref<const Vector3<scalar>>& eb1)
    {
        const Vector3<scalar> normal = (ea1 - ea0).cross(eb1 - eb0).normalized();
        return (eb0 - ea0).dot(normal) * normal;
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
        const Eigen::Ref<const Vector3<scalar>>& ea0,
        const Eigen::Ref<const Vector3<scalar>>& ea1,
        const Eigen::Ref<const Vector3<scalar>>& eb0,
        const Eigen::Ref<const Vector3<scalar>>& eb1,
        EdgeEdgeDistanceType dtype)
    {
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
            return point_line_closest_point_direction(eb0, ea0, ea1);

        case EdgeEdgeDistanceType::EA_EB1:
            return point_line_closest_point_direction(eb1, ea0, ea1);

        case EdgeEdgeDistanceType::EA0_EB:
            return -point_line_closest_point_direction(ea0, eb0, eb1);

        case EdgeEdgeDistanceType::EA1_EB:
            return -point_line_closest_point_direction(ea1, eb0, eb1);

        case EdgeEdgeDistanceType::EA_EB:
            return line_line_closest_point_direction<scalar>(ea0, ea1, eb0, eb1);

        default:
            throw std::invalid_argument(
                "Invalid distance type for edge-edge distance!");
        }
    }
}