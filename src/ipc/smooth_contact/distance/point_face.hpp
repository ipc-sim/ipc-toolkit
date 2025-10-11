#pragma once

#include "point_edge.hpp"

namespace ipc {

template <typename scalar>
scalar point_plane_sqr_distance(
    Eigen::ConstRef<Vector3<scalar>> p,
    Eigen::ConstRef<Vector3<scalar>> f0,
    Eigen::ConstRef<Vector3<scalar>> f1,
    Eigen::ConstRef<Vector3<scalar>> f2);

template <typename scalar>
scalar point_triangle_sqr_distance(
    Eigen::ConstRef<Vector3<scalar>> p,
    Eigen::ConstRef<Vector3<scalar>> t0,
    Eigen::ConstRef<Vector3<scalar>> t1,
    Eigen::ConstRef<Vector3<scalar>> t2,
    PointTriangleDistanceType dtype);

template <typename scalar>
Vector3<scalar> point_plane_closest_point_direction(
    Eigen::ConstRef<Vector3<scalar>> p,
    Eigen::ConstRef<Vector3<scalar>> f0,
    Eigen::ConstRef<Vector3<scalar>> f1,
    Eigen::ConstRef<Vector3<scalar>> f2);

std::tuple<Vector3d, Eigen::Matrix<double, 3, 12>>
point_plane_closest_point_direction_grad(
    Eigen::ConstRef<Vector3d> p,
    Eigen::ConstRef<Vector3d> t0,
    Eigen::ConstRef<Vector3d> t1,
    Eigen::ConstRef<Vector3d> t2);

std::tuple<Vector3d, Eigen::Matrix<double, 3, 12>, std::array<Matrix12d, 3>>
point_plane_closest_point_direction_hessian(
    Eigen::ConstRef<Vector3d> p,
    Eigen::ConstRef<Vector3d> t0,
    Eigen::ConstRef<Vector3d> t1,
    Eigen::ConstRef<Vector3d> t2);

template <typename scalar>
Vector3<scalar> point_triangle_closest_point_direction(
    Eigen::ConstRef<Vector3<scalar>> p,
    Eigen::ConstRef<Vector3<scalar>> t0,
    Eigen::ConstRef<Vector3<scalar>> t1,
    Eigen::ConstRef<Vector3<scalar>> t2,
    PointTriangleDistanceType dtype);

std::tuple<Vector3d, Eigen::Matrix<double, 3, 12>>
point_triangle_closest_point_direction_grad(
    Eigen::ConstRef<Vector3d> p,
    Eigen::ConstRef<Vector3d> t0,
    Eigen::ConstRef<Vector3d> t1,
    Eigen::ConstRef<Vector3d> t2,
    const PointTriangleDistanceType dtype = PointTriangleDistanceType::AUTO);

std::tuple<Vector3d, Eigen::Matrix<double, 3, 12>, std::array<Matrix12d, 3>>
point_triangle_closest_point_direction_hessian(
    Eigen::ConstRef<Vector3d> p,
    Eigen::ConstRef<Vector3d> t0,
    Eigen::ConstRef<Vector3d> t1,
    Eigen::ConstRef<Vector3d> t2,
    const PointTriangleDistanceType dtype = PointTriangleDistanceType::AUTO);
} // namespace ipc

#include "point_face.tpp"
