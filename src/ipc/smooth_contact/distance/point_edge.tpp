#pragma once

#include "point_edge.hpp"

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
            return (pos - L_ns(s) * t).squaredNorm();
        }
    }

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
            const scalar len = t.norm();
            t = t / len;

            const VectorMax3<scalar> pos = p - e0;
            const scalar s = pos.dot(t) / len;
            return pos - (L_ns(s) * len) * t;
        }
    }
}
