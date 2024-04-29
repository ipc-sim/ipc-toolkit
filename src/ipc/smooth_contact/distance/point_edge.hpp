#pragma once

#include <ipc/distance/distance_type.hpp>
#include <ipc/utils/math.hpp>
#include <ipc/smooth_contact/common.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/utils/AutodiffTypes.hpp>
#include <iostream>

namespace ipc {
template <typename scalar, int dim> class PointEdgeDistance {
public:
    static scalar point_point_sqr_distance(
        const Eigen::Ref<const Vector<scalar, dim>>& a,
        const Eigen::Ref<const Vector<scalar, dim>>& b);

    static scalar point_line_sqr_distance(
        const Eigen::Ref<const Vector<scalar, dim>>& p,
        const Eigen::Ref<const Vector<scalar, dim>>& e0,
        const Eigen::Ref<const Vector<scalar, dim>>& e1);

    static scalar point_edge_sqr_distance(
        const Eigen::Ref<const Vector<scalar, dim>>& p,
        const Eigen::Ref<const Vector<scalar, dim>>& e0,
        const Eigen::Ref<const Vector<scalar, dim>>& e1,
        const PointEdgeDistanceType dtype = PointEdgeDistanceType::AUTO);

    static Vector<scalar, dim> point_line_closest_point_direction(
        const Eigen::Ref<const Vector<scalar, dim>>& p,
        const Eigen::Ref<const Vector<scalar, dim>>& e0,
        const Eigen::Ref<const Vector<scalar, dim>>& e1);

    static Vector<scalar, dim> point_edge_closest_point_direction(
        const Eigen::Ref<const Vector<scalar, dim>>& p,
        const Eigen::Ref<const Vector<scalar, dim>>& e0,
        const Eigen::Ref<const Vector<scalar, dim>>& e1,
        const PointEdgeDistanceType dtype = PointEdgeDistanceType::AUTO);
};

template <int dim> class PointEdgeDistanceDerivatives {
public:
    static std::
        tuple<Vector<double, dim>, Eigen::Matrix<double, dim, dim * dim>>
        point_line_closest_point_direction_grad(
            const Eigen::Ref<const Vector<double, dim>>& p,
            const Eigen::Ref<const Vector<double, dim>>& e0,
            const Eigen::Ref<const Vector<double, dim>>& e1);

    static std::tuple<
        Vector<double, dim>,
        Eigen::Matrix<double, dim, dim * dim>,
        std::array<Eigen::Matrix<double, dim * dim, dim * dim>, dim>>
    point_line_closest_point_direction_hessian(
        const Eigen::Ref<const Vector<double, dim>>& p,
        const Eigen::Ref<const Vector<double, dim>>& e0,
        const Eigen::Ref<const Vector<double, dim>>& e1);

    static std::
        tuple<Vector<double, dim>, Eigen::Matrix<double, dim, dim * dim>>
        point_edge_closest_point_direction_grad(
            const Eigen::Ref<const Vector<double, dim>>& p,
            const Eigen::Ref<const Vector<double, dim>>& e0,
            const Eigen::Ref<const Vector<double, dim>>& e1,
            const PointEdgeDistanceType dtype = PointEdgeDistanceType::AUTO);

    static std::tuple<
        Vector<double, dim>,
        Eigen::Matrix<double, dim, dim * dim>,
        std::array<Eigen::Matrix<double, dim * dim, dim * dim>, dim>>
    point_edge_closest_point_direction_hessian(
        const Eigen::Ref<const Vector<double, dim>>& p,
        const Eigen::Ref<const Vector<double, dim>>& e0,
        const Eigen::Ref<const Vector<double, dim>>& e1,
        const PointEdgeDistanceType dtype = PointEdgeDistanceType::AUTO);
};
} // namespace ipc

#include "point_edge.tpp"
