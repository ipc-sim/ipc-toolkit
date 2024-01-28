#pragma once

#include <ipc/utils/distance_autodiff.hpp>

namespace ipc {

/// @brief 
/// @param v 
/// @param direc normalized, points to v
/// @param e0 
/// @param e1 
/// @param alpha 
/// @param beta 
/// @return 
template <class scalar>
inline scalar smooth_point2_term(
    const Eigen::Ref<const Vector2<scalar>>& v,
    const Eigen::Ref<const Vector2<scalar>>& direc,
    const Eigen::Ref<const Vector2<scalar>>& e0,
    const Eigen::Ref<const Vector2<scalar>>& e1,
    const double &alpha,
    const double &beta)
{
    const Vector2<scalar> t0 = (e0 - v).normalized(), t1 = (v - e1).normalized();

    const scalar tangent_term = smooth_heaviside<scalar>(direc.dot(t0), alpha, beta) *
                        smooth_heaviside<scalar>(-direc.dot(t1), alpha, beta);

    const scalar tmp = smooth_heaviside<scalar>(-cross2<scalar>(direc, t0), alpha, beta) + 
                         smooth_heaviside<scalar>(-cross2<scalar>(direc, t1), alpha, beta);
    const scalar normal_term = smooth_heaviside<scalar>(tmp - 1., alpha, 0);

    return tangent_term * normal_term * ((e0 - v).norm() + (e1 - v).norm()) / 2.;
}

inline bool smooth_point2_term_type(
    const Eigen::Ref<const Vector2<double>>& v,
    const Eigen::Ref<const Vector2<double>>& direc,
    const Eigen::Ref<const Vector2<double>>& e0,
    const Eigen::Ref<const Vector2<double>>& e1,
    const double &alpha,
    const double &beta)
{
    const Vector2<double> t0 = (e0 - v).normalized(), t1 = (v - e1).normalized();

    if (direc.dot(t0) <= -alpha || -direc.dot(t1) <= -alpha)
        return false;

    const double tmp = smooth_heaviside<double>(-cross2<double>(direc, t0), alpha, beta) + 
                         smooth_heaviside<double>(-cross2<double>(direc, t1), alpha, beta);
    if (tmp <= 1. - alpha)
        return false;

    return true;
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
    const double &alpha,
    const double &beta)
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
        const scalar tmp0 = smooth_heaviside<scalar>(direc.dot(t) / t.norm(), alpha, beta);

        tangent_term = tangent_term * tmp0;
        weight = weight + t.squaredNorm();

        normal_term = normal_term + smooth_heaviside<scalar>(-direc.dot(t_prev.cross(t).normalized()), alpha, beta);
        std::swap(t, t_prev);
    }

    // if constexpr (std::is_same<double, scalar>::value) 
    //     if (debug)
    //         logger().warn("normal {}, tangent {}", normal_term, tangent_term);

    return weight / 3 * tangent_term * smooth_heaviside<scalar>(normal_term - 1, alpha, 0);
}

inline bool smooth_point3_term_type(
    const Eigen::Ref<const RowVector3<double>>& v,
    const Eigen::Ref<const RowVector3<double>>& direc,
    const Eigen::Matrix<double, -1, 3> &neighbors,
    const double &alpha,
    const double &beta)
{
    RowVector3<double> t, t_prev;
    assert(neighbors.rows() > 2);

    bool tangent_term = true;
    double normal_term = 0;
    t_prev = neighbors.row(neighbors.rows()-1) - v;
    for (int a = 0; a < neighbors.rows(); a++)
    {
        t = neighbors.row(a) - v;
        tangent_term = tangent_term && ((direc.dot(t) / t.norm()) > - alpha);

        // logger().warn("normal term direction {}, should be negative", direc.dot(t_prev.cross(t).normalized()));

        normal_term += smooth_heaviside<double>(-direc.dot(t_prev.cross(t).normalized()), alpha, beta);
        std::swap(t, t_prev);
    }

    return  tangent_term && (normal_term > 1 - alpha);
}
}