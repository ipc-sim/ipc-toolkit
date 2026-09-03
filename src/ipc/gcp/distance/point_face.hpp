#pragma once

#include "point_edge.hpp"

namespace ipc {

template <typename T>
T point_plane_sqr_distance(
    Eigen::ConstRef<Eigen::Vector3<T>> p,
    Eigen::ConstRef<Eigen::Vector3<T>> f0,
    Eigen::ConstRef<Eigen::Vector3<T>> f1,
    Eigen::ConstRef<Eigen::Vector3<T>> f2)
{
    const Eigen::Vector3<T> normal = (f2 - f0).cross(f1 - f0);
    return Math<T>::sqr(normal.dot(p - f0)) / normal.squaredNorm();
}

template <typename scalar>
scalar point_triangle_sqr_distance(
    Eigen::ConstRef<Eigen::Vector3<scalar>> p,
    Eigen::ConstRef<Eigen::Vector3<scalar>> t0,
    Eigen::ConstRef<Eigen::Vector3<scalar>> t1,
    Eigen::ConstRef<Eigen::Vector3<scalar>> t2,
    PointTriangleDistanceType dtype)
{
    if (dtype == PointTriangleDistanceType::AUTO) {
        if constexpr (std::is_same<double, scalar>::value) {
            dtype = point_triangle_distance_type(p, t0, t1, t2);
        } else {
            Eigen::Vector3d p_, t0_, t1_, t2_;
            for (int d = 0; d < 3; d++) {
                p_(d) = p(d).val;
                t0_(d) = t0(d).val;
                t1_(d) = t1(d).val;
                t2_(d) = t2(d).val;
            }
            dtype = point_triangle_distance_type(p_, t0_, t1_, t2_);
        }
    }

    switch (dtype) {
    case PointTriangleDistanceType::P_T0: {
        return PointEdgeDistance<scalar, 3>::point_point_sqr_distance(p, t0);
    }

    case PointTriangleDistanceType::P_T1: {
        return PointEdgeDistance<scalar, 3>::point_point_sqr_distance(p, t1);
    }

    case PointTriangleDistanceType::P_T2: {
        return PointEdgeDistance<scalar, 3>::point_point_sqr_distance(p, t2);
    }

    case PointTriangleDistanceType::P_E0: {
        return PointEdgeDistance<scalar, 3>::point_line_sqr_distance(p, t0, t1);
    }

    case PointTriangleDistanceType::P_E1: {
        return PointEdgeDistance<scalar, 3>::point_line_sqr_distance(p, t1, t2);
    }

    case PointTriangleDistanceType::P_E2: {
        return PointEdgeDistance<scalar, 3>::point_line_sqr_distance(p, t2, t0);
    }

    case PointTriangleDistanceType::P_T: {
        return point_plane_sqr_distance<scalar>(p, t0, t1, t2);
    }

    default: {
        throw std::invalid_argument(
            "Invalid distance type for point-triangle distance!");
    }
    }
}

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
