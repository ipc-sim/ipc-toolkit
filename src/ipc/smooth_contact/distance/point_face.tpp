#pragma once

#include "point_face.hpp"

namespace ipc {
    template <typename scalar>
    scalar point_plane_sqr_distance(
        const Eigen::Ref<const Vector3<scalar>>& p,
        const Eigen::Ref<const Vector3<scalar>>& f0,
        const Eigen::Ref<const Vector3<scalar>>& f1,
        const Eigen::Ref<const Vector3<scalar>>& f2)
    {
        const Vector3<scalar> normal = (f2 - f0).cross(f1 - f0);
        return Math<scalar>::sqr(normal.dot(p - f0)) / normal.squaredNorm();
    }

    template <typename scalar>
    scalar point_triangle_sqr_distance(
        const Eigen::Ref<const Vector3<scalar>>& p,
        const Eigen::Ref<const Vector3<scalar>>& t0,
        const Eigen::Ref<const Vector3<scalar>>& t1,
        const Eigen::Ref<const Vector3<scalar>>& t2,
        PointTriangleDistanceType dtype)
    {
        if constexpr (std::is_same<double, scalar>::value)
            if (dtype == PointTriangleDistanceType::AUTO)
                dtype = point_triangle_distance_type(p, t0, t1, t2);

        switch (dtype) {
        case PointTriangleDistanceType::P_T0:
            return PointEdgeDistance<scalar, 3>::point_point_sqr_distance(p, t0);

        case PointTriangleDistanceType::P_T1:
            return PointEdgeDistance<scalar, 3>::point_point_sqr_distance(p, t1);

        case PointTriangleDistanceType::P_T2:
            return PointEdgeDistance<scalar, 3>::point_point_sqr_distance(p, t2);

        case PointTriangleDistanceType::P_E0:
            return PointEdgeDistance<scalar, 3>::point_line_sqr_distance(p, t0, t1);

        case PointTriangleDistanceType::P_E1:
            return PointEdgeDistance<scalar, 3>::point_line_sqr_distance(p, t1, t2);

        case PointTriangleDistanceType::P_E2:
            return PointEdgeDistance<scalar, 3>::point_line_sqr_distance(p, t2, t0);

        case PointTriangleDistanceType::P_T:
            return point_plane_sqr_distance<scalar>(p, t0, t1, t2);

        default:
            throw std::invalid_argument(
                "Invalid distance type for point-triangle distance!");
        }
    }

    template <typename scalar>
    Vector3<scalar> point_plane_closest_point_direction(
        const Eigen::Ref<const Vector3<scalar>>& p,
        const Eigen::Ref<const Vector3<scalar>>& f0,
        const Eigen::Ref<const Vector3<scalar>>& f1,
        const Eigen::Ref<const Vector3<scalar>>& f2)
    {
        const Vector3<scalar> normal = (f2 - f0).cross(f1 - f0);
        return (normal.dot(p - f0) / normal.squaredNorm()) * normal;
    }

    template <typename scalar>
    Vector3<scalar> point_triangle_closest_point_direction(
        const Eigen::Ref<const Vector3<scalar>>& p,
        const Eigen::Ref<const Vector3<scalar>>& t0,
        const Eigen::Ref<const Vector3<scalar>>& t1,
        const Eigen::Ref<const Vector3<scalar>>& t2,
        PointTriangleDistanceType dtype)
    {
        if constexpr (std::is_same<double, scalar>::value)
            if (dtype == PointTriangleDistanceType::AUTO)
                dtype = point_triangle_distance_type(p, t0, t1, t2);

        switch (dtype) {
        case PointTriangleDistanceType::P_T0:
            return p - t0;

        case PointTriangleDistanceType::P_T1:
            return p - t1;

        case PointTriangleDistanceType::P_T2:
            return p - t2;

        case PointTriangleDistanceType::P_E0:
            return PointEdgeDistance<scalar, 3>::point_line_closest_point_direction(p, t0, t1);

        case PointTriangleDistanceType::P_E1:
            return PointEdgeDistance<scalar, 3>::point_line_closest_point_direction(p, t1, t2);

        case PointTriangleDistanceType::P_E2:
            return PointEdgeDistance<scalar, 3>::point_line_closest_point_direction(p, t2, t0);

        case PointTriangleDistanceType::P_T:
            return point_plane_closest_point_direction<scalar>(p, t0, t1, t2);

        default:
            throw std::invalid_argument(
                "Invalid distance type for point-triangle distance!");
        }
    }
}