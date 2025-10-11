#pragma once

#include <ipc/distance/distance_type.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/smooth_contact/common.hpp>
#include <ipc/utils/autodiff_types.hpp>
#include <ipc/utils/math.hpp>

#include <iostream>

namespace ipc {
template <typename scalar, int dim> class PointEdgeDistance {
public:
    static scalar point_point_sqr_distance(
        Eigen::ConstRef<Vector<scalar, dim>> a,
        Eigen::ConstRef<Vector<scalar, dim>> b);

    static scalar point_line_sqr_distance(
        Eigen::ConstRef<Vector<scalar, dim>> p,
        Eigen::ConstRef<Vector<scalar, dim>> e0,
        Eigen::ConstRef<Vector<scalar, dim>> e1);

    static scalar point_edge_sqr_distance(
        Eigen::ConstRef<Vector<scalar, dim>> p,
        Eigen::ConstRef<Vector<scalar, dim>> e0,
        Eigen::ConstRef<Vector<scalar, dim>> e1,
        const PointEdgeDistanceType dtype = PointEdgeDistanceType::AUTO);

    static Vector<scalar, dim> point_line_closest_point_direction(
        Eigen::ConstRef<Vector<scalar, dim>> p,
        Eigen::ConstRef<Vector<scalar, dim>> e0,
        Eigen::ConstRef<Vector<scalar, dim>> e1);

    static Vector<scalar, dim> point_edge_closest_point_direction(
        Eigen::ConstRef<Vector<scalar, dim>> p,
        Eigen::ConstRef<Vector<scalar, dim>> e0,
        Eigen::ConstRef<Vector<scalar, dim>> e1,
        const PointEdgeDistanceType dtype = PointEdgeDistanceType::AUTO);
};

template <int dim> class PointEdgeDistanceDerivatives {
public:
    static std::
        tuple<Vector<double, dim>, Eigen::Matrix<double, dim, dim * dim>>
        point_line_closest_point_direction_grad(
            Eigen::ConstRef<Vector<double, dim>> p,
            Eigen::ConstRef<Vector<double, dim>> e0,
            Eigen::ConstRef<Vector<double, dim>> e1);

    static std::tuple<
        Vector<double, dim>,
        Eigen::Matrix<double, dim, dim * dim>,
        std::array<Eigen::Matrix<double, dim * dim, dim * dim>, dim>>
    point_line_closest_point_direction_hessian(
        Eigen::ConstRef<Vector<double, dim>> p,
        Eigen::ConstRef<Vector<double, dim>> e0,
        Eigen::ConstRef<Vector<double, dim>> e1);

    static std::
        tuple<Vector<double, dim>, Eigen::Matrix<double, dim, dim * dim>>
        point_edge_closest_point_direction_grad(
            Eigen::ConstRef<Vector<double, dim>> p,
            Eigen::ConstRef<Vector<double, dim>> e0,
            Eigen::ConstRef<Vector<double, dim>> e1,
            const PointEdgeDistanceType dtype = PointEdgeDistanceType::AUTO);

    static std::tuple<
        Vector<double, dim>,
        Eigen::Matrix<double, dim, dim * dim>,
        std::array<Eigen::Matrix<double, dim * dim, dim * dim>, dim>>
    point_edge_closest_point_direction_hessian(
        Eigen::ConstRef<Vector<double, dim>> p,
        Eigen::ConstRef<Vector<double, dim>> e0,
        Eigen::ConstRef<Vector<double, dim>> e1,
        const PointEdgeDistanceType dtype = PointEdgeDistanceType::AUTO);
};
} // namespace ipc

#include "point_edge.tpp"
