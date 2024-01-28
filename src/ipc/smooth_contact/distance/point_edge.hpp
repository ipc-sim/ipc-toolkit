#pragma once

#include <ipc/distance/distance_type.hpp>
#include <ipc/utils/math.hpp>

namespace ipc {
    template <typename scalar>
    scalar point_point_sqr_distance(
        const Eigen::Ref<const Vector3<scalar>>& a,
        const Eigen::Ref<const Vector3<scalar>>& b);

    template <typename scalar>
    scalar point_line_sqr_distance(
        const Eigen::Ref<const Vector3<scalar>>& p,
        const Eigen::Ref<const Vector3<scalar>>& e0,
        const Eigen::Ref<const Vector3<scalar>>& e1);

    template <typename scalar>
    scalar point_edge_sqr_distance(
        const Eigen::Ref<const VectorMax3<scalar>>& p,
        const Eigen::Ref<const VectorMax3<scalar>>& e0,
        const Eigen::Ref<const VectorMax3<scalar>>& e1,
        const PointEdgeDistanceType dtype = PointEdgeDistanceType::AUTO);

    template <typename scalar>
    VectorMax3<scalar> point_line_closest_point_direction(
        const Eigen::Ref<const VectorMax3<scalar>>& p,
        const Eigen::Ref<const VectorMax3<scalar>>& e0,
        const Eigen::Ref<const VectorMax3<scalar>>& e1);
    
    template <typename scalar>
    VectorMax3<scalar> point_edge_closest_point_direction(
        const Eigen::Ref<const VectorMax3<scalar>>& p,
        const Eigen::Ref<const VectorMax3<scalar>>& e0,
        const Eigen::Ref<const VectorMax3<scalar>>& e1,
        const PointEdgeDistanceType &dtype = PointEdgeDistanceType::AUTO);
}

#include "point_edge.tpp"
