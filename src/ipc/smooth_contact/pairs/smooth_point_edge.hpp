#pragma once

#include <ipc/smooth_contact/common.hpp>
#include <ipc/smooth_contact/distance/mollifier.hpp>
#include <ipc/smooth_contact/primitives/edge.hpp>
#include <ipc/smooth_contact/primitives/point2.hpp>

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
        direc = -direc / dist;

        const scalar mollifier_val = point_edge_mollifier<scalar>(p, e0, e1, dist*dist);
        return smooth_edge2_term<scalar>(direc, e1 - e0) * Math<scalar>::inv_barrier(dist / params.dhat, params.r) * mollifier_val * 
            smooth_point2_term<scalar>(p, direc, x0, x1, params.alpha, params.beta);
    }
}