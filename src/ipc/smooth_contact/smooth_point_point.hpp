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

/// @brief 
/// @tparam scalar 
/// @param direc normalized
/// @param v 
/// @param direc points to v
/// @param neighbors follow counter-clockwise order
/// @param params 
/// @return 
template <typename scalar>
inline scalar smooth_point3_term(
    const Eigen::Ref<const RowVector3<scalar>>& v,
    const Eigen::Ref<const RowVector3<scalar>>& direc,
    const Eigen::Matrix<scalar, -1, 3> &neighbors,
    const double &alpha)
{
    RowVector3<scalar> t, t_prev;
    assert(neighbors.rows() > 2);

    scalar tangent_term(1.);
    scalar weight(0.);
    scalar normal_term(0.);
    t_prev = neighbors.row(neighbors.rows()-1) - v;
    for (int a = 0; a < neighbors.rows(); a++)
    {
        t = neighbors.row(a) - v;
        tangent_term = tangent_term * smooth_heaviside<scalar>(direc.dot(t) / t.norm() / alpha);
        weight = weight + t.norm();

        normal_term = normal_term + smooth_heaviside<scalar>(-direc.dot(t_prev.cross(t).normalized()) / alpha);
        std::swap(t, t_prev);
    }

    return tangent_term * smooth_heaviside<scalar>(normal_term - 2) * weight / 3;
}

inline bool smooth_point3_term_type(
    const Eigen::Ref<const RowVector3<double>>& v,
    const Eigen::Ref<const RowVector3<double>>& direc,
    const Eigen::Matrix<double, -1, 3> &neighbors,
    const double &alpha)
{
    RowVector3<double> t, t_prev;
    assert(neighbors.rows() > 2);

    bool tangent_term = true;
    double normal_term = 0;
    t_prev = neighbors.row(neighbors.rows()-1) - v;
    for (int a = 0; a < neighbors.rows(); a++)
    {
        t = neighbors.row(a) - v;
        tangent_term = tangent_term && (direc.dot(t) / t.norm() / alpha > -1);

        normal_term += smooth_heaviside<double>(-direc.dot(t_prev.cross(t).normalized()) / alpha);
        std::swap(t, t_prev);
    }

    return tangent_term && (normal_term > 1);
}
// template <typename scalar>
// scalar smooth_point3_term(
//     const Eigen::Ref<const RowVector3<scalar>>& v,
//     const Eigen::Ref<const RowVector3<scalar>>& direc,
//     const Eigen::Matrix<scalar, -1, 3> &neighbors,
//     const double &alpha)
// {
//     RowVector3<scalar> t;
//     assert(neighbors.rows() > 2);

//     scalar tangent_term(1.);
//     scalar weight(0.);
//     for (int a = 0; a < neighbors.rows(); a++)
//     {
//         t = neighbors.row(a) - v;
//         tangent_term = tangent_term * smooth_heaviside<scalar>(direc.dot(t) / t.norm() / alpha);
//         weight = weight + t.squaredNorm();
//     }

//     return tangent_term * weight / neighbors.rows();
// }

// /// @brief 
// /// @param v 
// /// @param ray normalized
// /// @param neighbors counter-clockwise
// /// @return 
// bool is_outside_object(
//     const Eigen::Ref<const RowVector3<double>>& v,
//     const Eigen::Ref<const RowVector3<double>>& ray,
//     const Eigen::Matrix<double, -1, 3> &neighbors);

// bool smooth_point3_term_type(
//     const Eigen::Ref<const RowVector3<double>>& v,
//     const Eigen::Ref<const RowVector3<double>>& direc,
//     const Eigen::Matrix<double, -1, 3> &neighbors,
//     const double &alpha);

template <typename scalar>
scalar smooth_point_point_potential_3d(
    const Eigen::Ref<const RowVector3<scalar>>& va,
    const Eigen::Ref<const RowVector3<scalar>>& vb,
    const Eigen::Matrix<scalar, -1, 3> &ra,
    const Eigen::Matrix<scalar, -1, 3> &rb,
    const ParameterType &params)
{
    RowVector3<scalar> direc = va - vb;
    const scalar dist = direc.norm();
    direc = direc / dist;

    const scalar term_a = smooth_point3_term<scalar>(va, direc, ra, params.alpha);
    const scalar term_b = smooth_point3_term<scalar>(vb, -direc, rb, params.alpha);
    const scalar barrier = inv_barrier<scalar>(intpow(dist, 2) / params.eps, params.r);

    return term_a * term_b * barrier;
}
}