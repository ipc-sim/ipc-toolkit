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
scalar smooth_point_point_potential_2d(
    const Eigen::Ref<const Vector2<scalar>>& va,
    const Eigen::Ref<const Vector2<scalar>>& vb,
    const Eigen::Ref<const Vector2<scalar>>& ea0,
    const Eigen::Ref<const Vector2<scalar>>& ea1,
    const Eigen::Ref<const Vector2<scalar>>& eb0,
    const Eigen::Ref<const Vector2<scalar>>& eb1,
    const ParameterType &params)
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
    out = inv_barrier<scalar>(dist_sqr / intpow(params.eps,2), params.r) *
            smooth_heaviside<scalar>(-direc.dot(ta0) / dist / a0 / params.alpha) *
            smooth_heaviside<scalar>(direc.dot(ta1) / dist / a1 / params.alpha) *
            smooth_heaviside<scalar>(direc.dot(tb0) / dist / b0 / params.alpha) *
            smooth_heaviside<scalar>(direc.dot(-tb1) / dist / b1 / params.alpha) *
            (a0 + a1) * (b0 + b1) / 4.;

    // normal term
    scalar normal_term = (smooth_heaviside<scalar>(cross2<scalar>(direc, ta0) / dist / a0 / params.alpha) + 
                         smooth_heaviside<scalar>(cross2<scalar>(direc, ta1) / dist / a1 / params.alpha)) *
                        (smooth_heaviside<scalar>(-cross2<scalar>(direc, tb0) / dist / b0 / params.alpha) + 
                         smooth_heaviside<scalar>(-cross2<scalar>(direc, tb1) / dist / b1 / params.alpha));

    return out * normal_term;
}

template <typename scalar>
scalar smooth_point_point_potential_3d(
    const Eigen::Ref<const RowVector3<scalar>>& va,
    const Eigen::Ref<const RowVector3<scalar>>& vb,
    const Eigen::Matrix<scalar, -1, 3> &ra,
    const Eigen::Matrix<scalar, -1, 3> &rb,
    const ParameterType &params,
    const std::array<double, 2> &dhats)
{
    RowVector3<scalar> direc = va - vb;
    const scalar dist = direc.norm();
    direc = direc / dist;
    RowVector3<scalar> t, t_prev;

    assert(ra.rows() > 2);
    assert(rb.rows() > 2);

    scalar out_a = inv_barrier<scalar>(intpow(dist / dhats[0], 2), params.r);
    scalar normal_term_a(0.);
    scalar weight(0.);
    t_prev = ra.row(ra.rows()-1) - va;
    for (int a = 0; a < ra.rows(); a++)
    {
        t = ra.row(a) - va;
        out_a = out_a * smooth_heaviside<scalar>(direc.dot(t) / t.norm() / params.alpha);
        weight = weight + t.norm();

        normal_term_a = normal_term_a + smooth_heaviside<scalar>(-direc.dot(t_prev.cross(t).normalized()) / params.alpha);
        std::swap(t, t_prev);
    }
    out_a = out_a * weight / ra.rows();

    direc = -direc;
    scalar out_b = inv_barrier<scalar>(intpow(dist / dhats[1], 2), params.r);
    scalar normal_term_b(0.);
    weight = scalar(0.);
    t_prev = rb.row(rb.rows()-1) - vb;
    for (int b = 0; b < rb.rows(); b++)
    {
        t = rb.row(b) - vb;
        out_b = out_b * smooth_heaviside<scalar>(direc.dot(t) / t.norm() / params.alpha);
        weight = weight + t.norm();

        normal_term_b = normal_term_b + smooth_heaviside<scalar>(-direc.dot(t_prev.cross(t).normalized()) / params.alpha);
        std::swap(t, t_prev);
    }
    out_b = out_b * weight / rb.rows();

    return (out_a + out_b) * normal_term_a * normal_term_b;
}
}