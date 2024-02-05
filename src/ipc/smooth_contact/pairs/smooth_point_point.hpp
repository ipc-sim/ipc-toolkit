#pragma once

#include <ipc/smooth_contact/primitives/point2.hpp>
#include <ipc/smooth_contact/common.hpp>

namespace ipc {
/// @brief 
/// @tparam scalar 
/// @param va 
/// @param vb 
/// @param ea0 edge [va, ea0]
/// @param ea1 edge [ea1, va]
/// @param eb0 edge [vb, eb0]
/// @param eb1 edge [eb1, vb]
/// @param params 
/// @return 
template <typename scalar>
scalar smooth_point_point_potential_2d(
    const Eigen::Ref<const Vector2<scalar>>& va,
    const Eigen::Ref<const Vector2<scalar>>& vb,
    const Eigen::Ref<const Vector2<scalar>>& ea0,
    const Eigen::Ref<const Vector2<scalar>>& ea1,
    const Eigen::Ref<const Vector2<scalar>>& eb0,
    const Eigen::Ref<const Vector2<scalar>>& eb1,
    const ParameterType &params)
{
    Vector2<scalar> direc = vb - va;
    const scalar dist = direc.norm();
    direc = -direc / dist;

    return Math<scalar>::inv_barrier(dist / params.dhat, params.r) *
            smooth_point2_term<scalar>(va, -direc, ea0, ea1, params.alpha, params.beta) *
            smooth_point2_term<scalar>(vb,  direc, eb0, eb1, params.alpha, params.beta);
}

}