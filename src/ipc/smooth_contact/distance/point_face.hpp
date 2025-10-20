#pragma once

#include "point_edge.hpp"

namespace ipc {

template <typename T>
T point_plane_sqr_distance(
    Eigen::ConstRef<Eigen::Vector3<T>> p,
    Eigen::ConstRef<Eigen::Vector3<T>> f0,
    Eigen::ConstRef<Eigen::Vector3<T>> f1,
    Eigen::ConstRef<Eigen::Vector3<T>> f2);

template <typename T>
T point_triangle_sqr_distance(
    Eigen::ConstRef<Eigen::Vector3<T>> p,
    Eigen::ConstRef<Eigen::Vector3<T>> t0,
    Eigen::ConstRef<Eigen::Vector3<T>> t1,
    Eigen::ConstRef<Eigen::Vector3<T>> t2,
    PointTriangleDistanceType dtype);

template <typename T>
Eigen::Vector3<T> point_plane_closest_point_direction(
    Eigen::ConstRef<Eigen::Vector3<T>> p,
    Eigen::ConstRef<Eigen::Vector3<T>> f0,
    Eigen::ConstRef<Eigen::Vector3<T>> f1,
    Eigen::ConstRef<Eigen::Vector3<T>> f2);

std::tuple<Eigen::Vector3d, Eigen::Matrix<double, 3, 12>>
point_plane_closest_point_direction_grad(
    Eigen::ConstRef<Eigen::Vector3d> p,
    Eigen::ConstRef<Eigen::Vector3d> t0,
    Eigen::ConstRef<Eigen::Vector3d> t1,
    Eigen::ConstRef<Eigen::Vector3d> t2);

std::tuple<
    Eigen::Vector3d,
    Eigen::Matrix<double, 3, 12>,
    std::array<Matrix12d, 3>>
point_plane_closest_point_direction_hessian(
    Eigen::ConstRef<Eigen::Vector3d> p,
    Eigen::ConstRef<Eigen::Vector3d> t0,
    Eigen::ConstRef<Eigen::Vector3d> t1,
    Eigen::ConstRef<Eigen::Vector3d> t2);

template <typename T>
Eigen::Vector3<T> point_triangle_closest_point_direction(
    Eigen::ConstRef<Eigen::Vector3<T>> p,
    Eigen::ConstRef<Eigen::Vector3<T>> t0,
    Eigen::ConstRef<Eigen::Vector3<T>> t1,
    Eigen::ConstRef<Eigen::Vector3<T>> t2,
    PointTriangleDistanceType dtype);

std::tuple<Eigen::Vector3d, Eigen::Matrix<double, 3, 12>>
point_triangle_closest_point_direction_grad(
    Eigen::ConstRef<Eigen::Vector3d> p,
    Eigen::ConstRef<Eigen::Vector3d> t0,
    Eigen::ConstRef<Eigen::Vector3d> t1,
    Eigen::ConstRef<Eigen::Vector3d> t2,
    const PointTriangleDistanceType dtype = PointTriangleDistanceType::AUTO);

std::tuple<
    Eigen::Vector3d,
    Eigen::Matrix<double, 3, 12>,
    std::array<Matrix12d, 3>>
point_triangle_closest_point_direction_hessian(
    Eigen::ConstRef<Eigen::Vector3d> p,
    Eigen::ConstRef<Eigen::Vector3d> t0,
    Eigen::ConstRef<Eigen::Vector3d> t1,
    Eigen::ConstRef<Eigen::Vector3d> t2,
    const PointTriangleDistanceType dtype = PointTriangleDistanceType::AUTO);
} // namespace ipc

#include "point_face.tpp"
