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

std::tuple<Vector3d, Eigen::Matrix<double, 3, 12>>
point_plane_closest_point_direction_grad(
    const Eigen::Ref<const Vector3d>& p,
    const Eigen::Ref<const Vector3d>& t0,
    const Eigen::Ref<const Vector3d>& t1,
    const Eigen::Ref<const Vector3d>& t2);

std::tuple<Vector3d, Eigen::Matrix<double, 3, 12>, std::array<Matrix12d, 3>>
point_plane_closest_point_direction_hessian(
    const Eigen::Ref<const Vector3d>& p,
    const Eigen::Ref<const Vector3d>& t0,
    const Eigen::Ref<const Vector3d>& t1,
    const Eigen::Ref<const Vector3d>& t2);

template <typename scalar>
Vector3<scalar> point_triangle_closest_point_direction(
    const Eigen::Ref<const Vector3<scalar>>& p,
    const Eigen::Ref<const Vector3<scalar>>& t0,
    const Eigen::Ref<const Vector3<scalar>>& t1,
    const Eigen::Ref<const Vector3<scalar>>& t2,
    PointTriangleDistanceType dtype);

std::tuple<Vector3d, Eigen::Matrix<double, 3, 12>>
point_triangle_closest_point_direction_grad(
    const Eigen::Ref<const Vector3d>& p,
    const Eigen::Ref<const Vector3d>& t0,
    const Eigen::Ref<const Vector3d>& t1,
    const Eigen::Ref<const Vector3d>& t2,
    const PointTriangleDistanceType dtype = PointTriangleDistanceType::AUTO);

std::tuple<Vector3d, Eigen::Matrix<double, 3, 12>, std::array<Matrix12d, 3>>
point_triangle_closest_point_direction_hessian(
    const Eigen::Ref<const Vector3d>& p,
    const Eigen::Ref<const Vector3d>& t0,
    const Eigen::Ref<const Vector3d>& t1,
    const Eigen::Ref<const Vector3d>& t2,
    const PointTriangleDistanceType dtype = PointTriangleDistanceType::AUTO);
} // namespace ipc

#include "point_face.tpp"
