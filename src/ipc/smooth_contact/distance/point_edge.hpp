#pragma once

#include <ipc/distance/distance_type.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/math/math.hpp>
#include <ipc/smooth_contact/common.hpp>
#include <ipc/utils/autodiff_types.hpp>

#include <iostream>

namespace ipc {
template <typename T, int dim> class PointEdgeDistance {
public:
    using VectorNT = Eigen::Vector<T, dim>;

    PointEdgeDistance() = delete;
    PointEdgeDistance(const PointEdgeDistance&) = delete;
    PointEdgeDistance& operator=(const PointEdgeDistance&) = delete;

    static T point_point_sqr_distance(
        Eigen::ConstRef<VectorNT> a, Eigen::ConstRef<VectorNT> b);

    static T point_line_sqr_distance(
        Eigen::ConstRef<VectorNT> p,
        Eigen::ConstRef<VectorNT> e0,
        Eigen::ConstRef<VectorNT> e1);

    static T point_edge_sqr_distance(
        Eigen::ConstRef<VectorNT> p,
        Eigen::ConstRef<VectorNT> e0,
        Eigen::ConstRef<VectorNT> e1,
        const PointEdgeDistanceType dtype = PointEdgeDistanceType::AUTO);

    static VectorNT point_line_closest_point_direction(
        Eigen::ConstRef<VectorNT> p,
        Eigen::ConstRef<VectorNT> e0,
        Eigen::ConstRef<VectorNT> e1);

    static VectorNT point_edge_closest_point_direction(
        Eigen::ConstRef<VectorNT> p,
        Eigen::ConstRef<VectorNT> e0,
        Eigen::ConstRef<VectorNT> e1,
        const PointEdgeDistanceType dtype = PointEdgeDistanceType::AUTO);
};

template <int dim> class PointEdgeDistanceDerivatives {
public:
    using VectorNd = Eigen::Vector<double, dim>;
    using JacobianType =
        std::tuple<VectorNd, Eigen::Matrix<double, dim, 3 * dim>>;
    using HessianType = std::tuple<
        VectorNd,
        Eigen::Matrix<double, dim, 3 * dim>,
        std::array<Eigen::Matrix<double, 3 * dim, 3 * dim>, dim>>;

    PointEdgeDistanceDerivatives() = delete;
    PointEdgeDistanceDerivatives(const PointEdgeDistanceDerivatives&) = delete;
    PointEdgeDistanceDerivatives&
    operator=(const PointEdgeDistanceDerivatives&) = delete;

    static JacobianType point_line_closest_point_direction_grad(
        Eigen::ConstRef<VectorNd> p,
        Eigen::ConstRef<VectorNd> e0,
        Eigen::ConstRef<VectorNd> e1);

    static HessianType point_line_closest_point_direction_hessian(
        Eigen::ConstRef<VectorNd> p,
        Eigen::ConstRef<VectorNd> e0,
        Eigen::ConstRef<VectorNd> e1);

    static JacobianType point_edge_closest_point_direction_grad(
        Eigen::ConstRef<VectorNd> p,
        Eigen::ConstRef<VectorNd> e0,
        Eigen::ConstRef<VectorNd> e1,
        const PointEdgeDistanceType dtype = PointEdgeDistanceType::AUTO);

    static HessianType point_edge_closest_point_direction_hessian(
        Eigen::ConstRef<VectorNd> p,
        Eigen::ConstRef<VectorNd> e0,
        Eigen::ConstRef<VectorNd> e1,
        const PointEdgeDistanceType dtype = PointEdgeDistanceType::AUTO);
};
} // namespace ipc
