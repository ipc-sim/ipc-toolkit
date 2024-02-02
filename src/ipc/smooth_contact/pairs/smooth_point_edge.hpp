#pragma once

#include <ipc/smooth_contact/common.hpp>
#include <ipc/smooth_contact/distance/mollifier.hpp>
#include <ipc/smooth_contact/primitives/edge.hpp>
#include <ipc/smooth_contact/primitives/point.hpp>

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
    inline scalar smooth_point_edge_potential_single_point(
        const Eigen::Ref<const Vector2<scalar>>& p,
        const Eigen::Ref<const Vector2<scalar>>& e0,
        const Eigen::Ref<const Vector2<scalar>>& e1,
        const Eigen::Ref<const Vector2<scalar>>& x0,
        const Eigen::Ref<const Vector2<scalar>>& x1,
        const ParameterType &params)
    {
        Vector2<scalar> direc = PointEdgeDistance<scalar, 2>::point_edge_closest_point_direction(p, e0, e1, PointEdgeDistanceType::AUTO);
        const scalar dist = direc.norm();
        direc = direc / dist;

        const scalar mollifier_val = point_edge_mollifier<scalar>(p, e0, e1, dist*dist);
        return smooth_edge2_term<scalar>(direc, e1 - e0) * Math<scalar>::inv_barrier(dist / params.dhat, params.r) * mollifier_val * 
            smooth_point2_term<scalar>(p, direc, x0, x1, params.alpha, params.beta);
    }

    inline bool smooth_point_edge_potential_single_point_3d_type(
        const Eigen::Ref<const Vector3<double>>& p,
        const Eigen::Matrix<double, -1, 3> &neighbors,
        const Eigen::Ref<const Vector3<double>>& e0,
        const Eigen::Ref<const Vector3<double>>& e1,
        const Eigen::Ref<const Vector3<double>>& f0,
        const Eigen::Ref<const Vector3<double>>& f1,
        const ParameterType &params,
        ORIENTATION_TYPES &point_otypes,
        ORIENTATION_TYPES &edge_otypes)
    {
        Vector3<double> direc = PointEdgeDistance<double, 3>::point_edge_closest_point_direction(p, e0, e1, PointEdgeDistanceType::AUTO); // from edge a to edge b
        const double dist = direc.norm();
        if (dist >= params.dhat)
            return false;
        
        direc = direc / dist;

        const bool edge_term = smooth_edge3_term_type(direc, e0, e1, f0, f1, params.alpha, params.beta, edge_otypes);
        const bool vert_term = smooth_point3_term_type(p, direc, neighbors, params.alpha, params.beta, point_otypes);

        return edge_term && vert_term;
    }

    template <typename scalar>
    scalar smooth_point_edge_potential_single_point_3d(
        const Eigen::Ref<const Vector3<scalar>>& p,
        const Eigen::Ref<const Eigen::Matrix<scalar, -1, 3>>& neighbors,
        const Eigen::Ref<const Vector3<scalar>>& e0,
        const Eigen::Ref<const Vector3<scalar>>& e1,
        const Eigen::Ref<const Vector3<scalar>>& f0,
        const Eigen::Ref<const Vector3<scalar>>& f1,
        const ParameterType &params,
        const ORIENTATION_TYPES &point_otypes,
        const ORIENTATION_TYPES &edge_otypes)
    {
        Vector3<scalar> direc = PointEdgeDistance<scalar, 3>::point_edge_closest_point_direction(p, e0, e1, PointEdgeDistanceType::AUTO); // from edge a to edge b
        const scalar dist = direc.norm();
        direc = direc / dist;

        const scalar edge_term = smooth_edge3_term<scalar>(direc, e0, e1, f0, f1, params.alpha, params.beta, edge_otypes);
        const scalar barrier = Math<scalar>::inv_barrier(dist / params.dhat, params.r);
        const scalar mollifier_val = point_edge_mollifier<scalar>(p, e0, e1, dist*dist);
        const scalar vert_term = smooth_point3_term<scalar>(p, direc, neighbors, params.alpha, params.beta, point_otypes);

        return edge_term * mollifier_val * vert_term * barrier;
    }
}