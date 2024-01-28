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
    const double &beta,
    const ORIENTATION_TYPES &otypes)
{
    RowVector3<scalar> t, t_prev;
    assert(neighbors.rows() > 2);
    assert(otypes.size() == neighbors.rows());

    scalar tangent_term(1.);
    scalar weight(0.);
    scalar normal_term(0.);
    t_prev = neighbors.row(neighbors.rows()-1) - v;
    for (int a = 0; a < neighbors.rows(); a++)
    {
        t = neighbors.row(a) - v;
        if (otypes.tangent_type(a) == HEAVISIDE_TYPE::VARIANT)
        {
            scalar tmp0 = smooth_heaviside<scalar>(direc.dot(t) / t.norm(), alpha, beta);
            tangent_term = tangent_term * tmp0;
        }
        else if (otypes.tangent_type(a) == HEAVISIDE_TYPE::ZERO)
        {
            tangent_term = scalar(0.);
            break;
        }

        // if normal_term >= 1, the term depending on it becomes constant 1
        if (normal_term < 1 && otypes.normal_type(a) != HEAVISIDE_TYPE::ZERO)
            normal_term = normal_term + (otypes.normal_type(a) == HEAVISIDE_TYPE::ONE ? scalar(1.) : smooth_heaviside<scalar>(-direc.dot(t_prev.cross(t).normalized()), alpha, beta));
        
        weight = weight + t.squaredNorm();
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
    const double &beta,
    ORIENTATION_TYPES &otypes)
{
    RowVector3<double> t, t_prev;
    assert(neighbors.rows() > 2);
    otypes.set_size(neighbors.rows());

    bool tangent_term = true;
    double normal_term = 0;
    t_prev = neighbors.row(neighbors.rows()-1) - v;
    for (int a = 0; a < neighbors.rows(); a++)
    {
        t = neighbors.row(a) - v;
        otypes.tangent_type(a) = otypes.compute_type(direc.dot(t) / t.norm(), alpha, beta);
        if (otypes.tangent_type(a) == HEAVISIDE_TYPE::ZERO)
            return false;

        const double tmp = -direc.dot(t_prev.cross(t).normalized());
        otypes.normal_type(a) = otypes.compute_type(tmp, alpha, beta);
        normal_term += smooth_heaviside<double>(tmp, alpha, beta);

        std::swap(t, t_prev);
    }

    return  tangent_term && (normal_term > 1 - alpha);
}
}