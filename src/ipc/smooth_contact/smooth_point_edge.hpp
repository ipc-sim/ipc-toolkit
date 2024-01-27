#pragma once

#include <ipc/distance/distance_type.hpp>
#include <ipc/collisions/collision.hpp>
#include "smooth_edge_edge.hpp"
#include "smooth_point_point.hpp"

namespace ipc {

    /// @brief 
    /// @tparam scalar 
    /// @param p 
    /// @param e0 
    /// @param e1 edge [e0, e1]
    /// @param t0 edge [p, t0]
    /// @param t1 edge [t1, p]
    /// @param params 
    /// @return 
    template <typename scalar>
    scalar smooth_point_edge_potential_single_point(
        const Eigen::Ref<const Vector2<scalar>>& p,
        const Eigen::Ref<const Vector2<scalar>>& e0,
        const Eigen::Ref<const Vector2<scalar>>& e1,
        const Eigen::Ref<const Vector2<scalar>>& x0,
        const Eigen::Ref<const Vector2<scalar>>& x1,
        const ParameterType &params)
    {
        Vector2<scalar> tangent = e1 - e0;
        const scalar len = tangent.norm();
        tangent = tangent / len;

        Vector2<scalar> direc = point_edge_closest_point_direction<scalar>(p, e0, e1, PointEdgeDistanceType::AUTO);
        const scalar dist_sqr = direc.squaredNorm();
        direc = direc / sqrt(dist_sqr);
        const scalar Phi = 1 - cross2<scalar>(direc, tangent);
        const scalar mollifier_val = edge_mollifier<scalar>(p, e0, e1, dist_sqr);

        return len * cubic_spline(Phi / params.alpha) * inv_barrier(dist_sqr / params.eps, params.r) * 
                mollifier_val * smooth_point2_term<scalar>(p, direc, x0, x1, params.alpha, params.beta);
    }

    inline bool smooth_point_edge_potential_single_point_3d_type(
        const Eigen::Ref<const Vector3<double>>& p,
        const Eigen::Matrix<double, -1, 3> &neighbors,
        const Eigen::Ref<const Vector3<double>>& e0,
        const Eigen::Ref<const Vector3<double>>& e1,
        const Eigen::Ref<const Vector3<double>>& f0,
        const Eigen::Ref<const Vector3<double>>& f1,
        const ParameterType &params)
    {
        Vector3<double> direc = point_edge_closest_point_direction<double>(p, e0, e1, PointEdgeDistanceType::AUTO); // from edge a to edge b
        const double dist_sqr = direc.squaredNorm();
        if (dist_sqr >= params.eps)
            return false;
        
        direc = direc / sqrt(dist_sqr);

        const bool edge_term = smooth_edge3_term_type(direc, e0, e1, f0, f1, params.alpha, params.beta);
        const bool vert_term = smooth_point3_term_type(p, direc, neighbors, params.alpha, params.beta);

        return edge_term && vert_term;
    }

    template <typename scalar>
    scalar smooth_point_edge_potential_single_point_3d(
        const Eigen::Ref<const Vector3<scalar>>& p,
        const Eigen::Matrix<scalar, -1, 3> &neighbors,
        const Eigen::Ref<const Vector3<scalar>>& e0,
        const Eigen::Ref<const Vector3<scalar>>& e1,
        const Eigen::Ref<const Vector3<scalar>>& f0,
        const Eigen::Ref<const Vector3<scalar>>& f1,
        const ParameterType &params)
    {
        Vector3<scalar> direc = point_edge_closest_point_direction<scalar>(p, e0, e1, PointEdgeDistanceType::AUTO); // from edge a to edge b
        const scalar dist_sqr = direc.squaredNorm();
        direc = direc / sqrt(dist_sqr);

        const scalar edge_term = smooth_edge3_term<scalar>(direc, e0, e1, f0, f1, params.alpha, params.beta);
        const scalar barrier = inv_barrier<scalar>(dist_sqr / params.eps, params.r);
        const scalar mollifier_val = edge_mollifier<scalar>(p, e0, e1, dist_sqr);
        const scalar vert_term = smooth_point3_term<scalar>(p, direc, neighbors, params.alpha, params.beta);

        return edge_term * mollifier_val * vert_term * barrier;
    }
}