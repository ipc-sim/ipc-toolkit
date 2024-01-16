#pragma once

#include <ipc/distance/distance_type.hpp>
#include <ipc/collisions/collision.hpp>
#include <ipc/utils/math.hpp>

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
scalar smooth_point_point_potential(
    const Eigen::Ref<const Vector2<scalar>>& va,
    const Eigen::Ref<const Vector2<scalar>>& vb,
    const Eigen::Ref<const Vector2<scalar>>& ea0,
    const Eigen::Ref<const Vector2<scalar>>& ea1,
    const Eigen::Ref<const Vector2<scalar>>& eb0,
    const Eigen::Ref<const Vector2<scalar>>& eb1,
    const ParameterType &params,
    const std::array<double, 2> &dhats)
{
    const Vector2<scalar> direc = vb - va;
    const scalar dist_sqr = direc.squaredNorm();
    const scalar dist = sqrt(dist_sqr);

    const Vector2<scalar> ta0 = ea0 - va, ta1 = va - ea1;
    const Vector2<scalar> tb0 = eb0 - vb, tb1 = vb - eb1;

    const scalar a0 = ta0.norm(), a1 = ta1.norm();
    const scalar b0 = tb0.norm(), b1 = tb1.norm();

    scalar out(0.);
    // tangent term
    out += inv_barrier<scalar>(dist_sqr / intpow(dhats[0],2), params.r) *
            smooth_heaviside<scalar>(-direc.dot(ta0) / dist / a0 / params.alpha) *
            smooth_heaviside<scalar>(direc.dot(ta1) / dist / a1 / params.alpha) *
            (a0 + a1) / 2.;
    out += inv_barrier<scalar>(dist_sqr / intpow(dhats[1],2), params.r) *
            smooth_heaviside<scalar>(direc.dot(tb0) / dist / b0 / params.alpha) *
            smooth_heaviside<scalar>(direc.dot(-tb1) / dist / b1 / params.alpha) *
            (b0 + b1) / 2.;

    // normal term
    scalar normal_term = (smooth_heaviside<scalar>(cross2<scalar>(direc, ta0) / dist / a0 / params.alpha) + 
                         smooth_heaviside<scalar>(cross2<scalar>(direc, ta1) / dist / a1 / params.alpha)) *
                        (smooth_heaviside<scalar>(-cross2<scalar>(direc, tb0) / dist / b0 / params.alpha) + 
                         smooth_heaviside<scalar>(-cross2<scalar>(direc, tb1) / dist / b1 / params.alpha));

    return out * normal_term;
}
}