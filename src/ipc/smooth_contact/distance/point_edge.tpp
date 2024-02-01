#pragma once

#include "point_edge.hpp"
#include <ipc/distance/point_edge.hpp>
#include <ipc/utils/AutodiffTypes.hpp>
#include <iostream>

namespace ipc {
    template <typename scalar>
    scalar point_point_sqr_distance(
        const Eigen::Ref<const Vector3<scalar>>& a,
        const Eigen::Ref<const Vector3<scalar>>& b)
    {
        return (a - b).squaredNorm();
    }

    template <typename scalar>
    scalar point_line_sqr_distance(
        const Eigen::Ref<const Vector3<scalar>>& p,
        const Eigen::Ref<const Vector3<scalar>>& e0,
        const Eigen::Ref<const Vector3<scalar>>& e1)
    {
        return (e0 - p).cross(e1 - p).squaredNorm() / (e1 - e0).squaredNorm();
    }

    template <typename scalar>
    scalar point_edge_sqr_distance(
        const Eigen::Ref<const VectorMax3<scalar>>& p,
        const Eigen::Ref<const VectorMax3<scalar>>& e0,
        const Eigen::Ref<const VectorMax3<scalar>>& e1,
        const PointEdgeDistanceType dtype)
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

    // template <int dim, int max_dim>
    // ADHessian<dim, max_dim> point_edge_sqr_distance(
    //     const Eigen::Ref<const VectorMax3<ADHessian<dim, max_dim>>>& p,
    //     const Eigen::Ref<const VectorMax3<ADHessian<dim, max_dim>>>& e0,
    //     const Eigen::Ref<const VectorMax3<ADHessian<dim, max_dim>>>& e1,
    //     const PointEdgeDistanceType dtype)
    // {
    //     const int size = p.size();
    //     VectorMax9d vec(size * 3);
    //     vec.head(size) = AutoDiffToDouble(p);
    //     vec.segment(size, size) = AutoDiffToDouble(e0);
    //     vec.tail(size) = AutoDiffToDouble(e1);

    //     double val = point_edge_distance(vec.head(size), vec.segment(size, size), vec.tail(size), dtype);
    //     auto grad = point_edge_distance_gradient(vec.head(size), vec.segment(size, size), vec.tail(size), dtype);
    //     auto hess = point_edge_distance_hessian(vec.head(size), vec.segment(size, size), vec.tail(size), dtype);

    //     ADHessian<dim, max_dim> out;
    //     out.value = val;

    //     MatrixMax9d G(size * 3, size * 3);
    //     for (int d = 0; d < size; d++)
    //     {
    //         G.col(d) = p(d).grad;
    //         G.col(size+d) = e0(d).grad;
    //         G.col(2*size+d) = e1(d).grad;
    //     }
    //     out.grad = G * grad;

    //     for (int d = 0; d < size; d++)
    //     {
    //         out.hess += p(d).hess * grad(d);
    //         out.hess += e0(d).hess * grad(size+d);
    //         out.hess += e1(d).hess * grad(2*size+d);
    //     }
    //     out.hess += G * hess * G.transpose();
    //     std::cout << "point-edge dist hess fast\n";

    //     return out;
    // }

    template <typename scalar>
    VectorMax3<scalar> point_line_closest_point_direction(
        const Eigen::Ref<const VectorMax3<scalar>>& p,
        const Eigen::Ref<const VectorMax3<scalar>>& e0,
        const Eigen::Ref<const VectorMax3<scalar>>& e1)
    {
        const VectorMax3<scalar> d = p - e0;
        const VectorMax3<scalar> t = e1 - e0;
        return d - (d.dot(t) / t.squaredNorm()) * t;
    }

    template <typename scalar>
    VectorMax3<scalar> point_edge_closest_point_direction(
        const Eigen::Ref<const VectorMax3<scalar>>& p,
        const Eigen::Ref<const VectorMax3<scalar>>& e0,
        const Eigen::Ref<const VectorMax3<scalar>>& e1,
        const PointEdgeDistanceType &dtype)
    {
        switch (dtype)
        {
        case PointEdgeDistanceType::P_E:
            return point_line_closest_point_direction<scalar>(p, e0, e1);
        case PointEdgeDistanceType::P_E0:
            return p - e0;
        case PointEdgeDistanceType::P_E1:
            return p - e1;
        case PointEdgeDistanceType::AUTO:
        default:
            VectorMax3<scalar> t = e1 - e0;
            const VectorMax3<scalar> pos = p - e0;
            const scalar s = pos.dot(t) / t.squaredNorm();
            return pos - Math<scalar>::L_ns(s) * t;
        }
    }
}
