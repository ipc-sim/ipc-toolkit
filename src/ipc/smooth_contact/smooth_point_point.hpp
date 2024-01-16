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

    const Vector2<scalar> ta0 = ea0 - va, ta1 = va - ea1;
    const Vector2<scalar> tb0 = eb0 - vb, tb1 = vb - eb1;

    scalar out(0.);
    out += inv_barrier<scalar>(dist_sqr / dhats[0], params.r) *
            smooth_heaviside<scalar>(-direc.dot(ta0) / sqrt(dist_sqr * ta0.squaredNorm()) / params.alpha) *
            smooth_heaviside<scalar>(-direc.dot(ta1) / sqrt(dist_sqr * ta1.squaredNorm()) / params.alpha) *
            (ta0.norm() + ta1.norm()) / 2.;
    out += inv_barrier<scalar>(dist_sqr / dhats[1], params.r) *
            smooth_heaviside<scalar>(direc.dot(tb0) / sqrt(dist_sqr * tb0.squaredNorm()) / params.alpha) *
            smooth_heaviside<scalar>(direc.dot(tb1) / sqrt(dist_sqr * tb1.squaredNorm()) / params.alpha) *
            (tb0.norm() + tb1.norm()) / 2.;
    return out;
}
}