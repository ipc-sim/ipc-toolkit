#pragma once

#include "point_edge.hpp"

namespace ipc {

    template <typename scalar>
    scalar point_plane_sqr_distance(
        const Eigen::Ref<const Vector3<scalar>>& p,
        const Eigen::Ref<const Vector3<scalar>>& f0,
        const Eigen::Ref<const Vector3<scalar>>& f1,
        const Eigen::Ref<const Vector3<scalar>>& f2);

    template <typename scalar>
    scalar point_triangle_sqr_distance(
        const Eigen::Ref<const Vector3<scalar>>& p,
        const Eigen::Ref<const Vector3<scalar>>& t0,
        const Eigen::Ref<const Vector3<scalar>>& t1,
        const Eigen::Ref<const Vector3<scalar>>& t2,
        PointTriangleDistanceType dtype);

    template <typename scalar>
    Vector3<scalar> point_plane_closest_point_direction(
        const Eigen::Ref<const Vector3<scalar>>& p,
        const Eigen::Ref<const Vector3<scalar>>& f0,
        const Eigen::Ref<const Vector3<scalar>>& f1,
        const Eigen::Ref<const Vector3<scalar>>& f2);

    template <typename scalar>
    Vector3<scalar> point_triangle_closest_point_direction(
        const Eigen::Ref<const Vector3<scalar>>& p,
        const Eigen::Ref<const Vector3<scalar>>& t0,
        const Eigen::Ref<const Vector3<scalar>>& t1,
        const Eigen::Ref<const Vector3<scalar>>& t2,
        PointTriangleDistanceType dtype);
}

#include "point_face.tpp"
