#pragma once

#include <ipc/distance/distance_type.hpp>
#include <ipc/utils/math.hpp>

#include <ipc/distance/point_edge.hpp>
#include <ipc/utils/AutodiffTypes.hpp>
#include <iostream>

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

    // template <typename scalar>
    // scalar point_edge_sqr_distance(
    //     const Eigen::Ref<const VectorMax3<scalar>>& p,
    //     const Eigen::Ref<const VectorMax3<scalar>>& e0,
    //     const Eigen::Ref<const VectorMax3<scalar>>& e1,
    //     const PointEdgeDistanceType dtype = PointEdgeDistanceType::AUTO)
    // {
    //     switch (dtype)
    //     {
    //     case PointEdgeDistanceType::P_E:
    //         return point_line_sqr_distance<scalar>(p, e0, e1);
    //     case PointEdgeDistanceType::P_E0:
    //         return point_point_sqr_distance<scalar>(p, e0);
    //     case PointEdgeDistanceType::P_E1:
    //         return point_point_sqr_distance<scalar>(p, e1);
    //     case PointEdgeDistanceType::AUTO:
    //     default:
    //         const VectorMax3<scalar> t = e1 - e0;
    //         const VectorMax3<scalar> pos = p - e0;
    //         const scalar s = pos.dot(t) / t.squaredNorm();
    //         return (pos - Math<scalar>::L_ns(s) * t).squaredNorm();
    //     }
    // }

    template <typename scalar>
    scalar point_edge_sqr_distance(
        const Eigen::Ref<const VectorMax3<scalar>>& p,
        const Eigen::Ref<const VectorMax3<scalar>>& e0,
        const Eigen::Ref<const VectorMax3<scalar>>& e1,
        const PointEdgeDistanceType dtype = PointEdgeDistanceType::AUTO)
    {
        if constexpr (std::is_same<scalar, ADHessian<12>>::value)
        {
            auto pd = autodiff_to_double<scalar, 3, 1, 3>(p);
            auto e0d = autodiff_to_double<scalar, 3, 1, 3>(e0);
            auto e1d = autodiff_to_double<scalar, 3, 1, 3>(e1);
            const double val = point_edge_distance(pd, e0d, e1d, dtype);
            const Vector9d grad = point_edge_distance_gradient(pd, e0d, e1d, dtype);
            const Matrix9d hess = point_edge_distance_hessian(pd, e0d, e1d, dtype);
            Eigen::Matrix<double, 12, 9> G;
            G << p(0).getGradient(), p(1).getGradient(), p(2).getGradient(),
                e0(0).getGradient(), e0(1).getGradient(), e0(2).getGradient(),
                e1(0).getGradient(), e1(1).getGradient(), e1(2).getGradient();
            Eigen::Matrix<double, 12, 12> H = G * hess * G.transpose();
            for (int i = 0; i < 3; i++)
            {
                H += p(i).getHessian() * grad(i);
                H += e0(i).getHessian() * grad(3+i);
                H += e1(i).getHessian() * grad(6+i);
            }
            return scalar(val, G * grad, H);
        }
        else
        {
            switch (dtype)
            {
            case PointEdgeDistanceType::P_E:
                return point_line_sqr_distance<scalar>(p, e0, e1);
            case PointEdgeDistanceType::P_E0:
                return point_point_sqr_distance<scalar>(p, e0);
            case PointEdgeDistanceType::P_E1:
                return point_point_sqr_distance<scalar>(p, e1);
            case PointEdgeDistanceType::AUTO:
            default:
                const VectorMax3<scalar> t = e1 - e0;
                const VectorMax3<scalar> pos = p - e0;
                const scalar s = pos.dot(t) / t.squaredNorm();
                return (pos - Math<scalar>::L_ns(s) * t).squaredNorm();
            }
        }
    }

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
