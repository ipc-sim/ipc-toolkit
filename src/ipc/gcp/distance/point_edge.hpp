#pragma once

#include <ipc/distance/distance_type.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/distance/point_line.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/math/math.hpp>
#include <ipc/gcp/common.hpp>
#include <ipc/utils/autodiff_types.hpp>

namespace ipc {
template <typename T, int dim> class PointEdgeDistance {
public:
    using VectorNT = Eigen::Vector<T, dim>;

    PointEdgeDistance() = delete;
    PointEdgeDistance(const PointEdgeDistance&) = delete;
    PointEdgeDistance& operator=(const PointEdgeDistance&) = delete;

    static T point_point_sqr_distance(
        Eigen::ConstRef<Eigen::Vector<T, dim>> a,
        Eigen::ConstRef<Eigen::Vector<T, dim>> b)
    {
        return (a - b).squaredNorm();
    }

    static T point_line_sqr_distance(
        Eigen::ConstRef<Eigen::Vector<T, dim>> p,
        Eigen::ConstRef<Eigen::Vector<T, dim>> e0,
        Eigen::ConstRef<Eigen::Vector<T, dim>> e1)
    {
        if constexpr (dim == 2) {
            return Math<T>::sqr(Math<T>::cross2(e0 - p, e1 - p))
                / (e1 - e0).squaredNorm();
        } else {
            return (e0 - p).cross(e1 - p).squaredNorm()
                / (e1 - e0).squaredNorm();
        }
    }

    static T point_edge_sqr_distance(
        Eigen::ConstRef<Eigen::Vector<T, dim>> p,
        Eigen::ConstRef<Eigen::Vector<T, dim>> e0,
        Eigen::ConstRef<Eigen::Vector<T, dim>> e1,
        const PointEdgeDistanceType dtype = PointEdgeDistanceType::AUTO)
    {
        switch (dtype) {
        case PointEdgeDistanceType::P_E:
            return point_line_sqr_distance(p, e0, e1);
        case PointEdgeDistanceType::P_E0:
            return point_point_sqr_distance(p, e0);
        case PointEdgeDistanceType::P_E1:
            return point_point_sqr_distance(p, e1);
        case PointEdgeDistanceType::AUTO:
        default:
            const Eigen::Vector<T, dim> t = e1 - e0;
            const Eigen::Vector<T, dim> pos = p - e0;
            const T s = pos.dot(t) / t.squaredNorm();
            return (pos - Math<T>::l_ns(s) * t).squaredNorm();
        }
    }

    static Eigen::Vector<T, dim> point_line_closest_point_direction(
        Eigen::ConstRef<Eigen::Vector<T, dim>> p,
        Eigen::ConstRef<Eigen::Vector<T, dim>> e0,
        Eigen::ConstRef<Eigen::Vector<T, dim>> e1);

    static Eigen::Vector<T, dim> point_edge_closest_point_direction(
        Eigen::ConstRef<Eigen::Vector<T, dim>> p,
        Eigen::ConstRef<Eigen::Vector<T, dim>> e0,
        Eigen::ConstRef<Eigen::Vector<T, dim>> e1,
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

    static std::
        tuple<Eigen::Vector<double, dim>, Eigen::Matrix<double, dim, 3 * dim>>
        point_line_closest_point_direction_grad(
            Eigen::ConstRef<Eigen::Vector<double, dim>> p,
            Eigen::ConstRef<Eigen::Vector<double, dim>> e0,
            Eigen::ConstRef<Eigen::Vector<double, dim>> e1);

    static std::tuple<
        Eigen::Vector<double, dim>,
        Eigen::Matrix<double, dim, 3 * dim>,
        std::array<Eigen::Matrix<double, 3 * dim, 3 * dim>, dim>>
    point_line_closest_point_direction_hessian(
        Eigen::ConstRef<Eigen::Vector<double, dim>> p,
        Eigen::ConstRef<Eigen::Vector<double, dim>> e0,
        Eigen::ConstRef<Eigen::Vector<double, dim>> e1);

    static std::
        tuple<Eigen::Vector<double, dim>, Eigen::Matrix<double, dim, 3 * dim>>
        point_edge_closest_point_direction_grad(
            Eigen::ConstRef<Eigen::Vector<double, dim>> p,
            Eigen::ConstRef<Eigen::Vector<double, dim>> e0,
            Eigen::ConstRef<Eigen::Vector<double, dim>> e1,
            const PointEdgeDistanceType dtype = PointEdgeDistanceType::AUTO);

    static std::tuple<
        Eigen::Vector<double, dim>,
        Eigen::Matrix<double, dim, 3 * dim>,
        std::array<Eigen::Matrix<double, 3 * dim, 3 * dim>, dim>>
    point_edge_closest_point_direction_hessian(
        Eigen::ConstRef<Eigen::Vector<double, dim>> p,
        Eigen::ConstRef<Eigen::Vector<double, dim>> e0,
        Eigen::ConstRef<Eigen::Vector<double, dim>> e1,
        const PointEdgeDistanceType dtype = PointEdgeDistanceType::AUTO);
};
} // namespace ipc
