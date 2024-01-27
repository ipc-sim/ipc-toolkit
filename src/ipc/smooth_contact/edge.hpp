#pragma once

#include <ipc/utils/distance_autodiff.hpp>

namespace ipc {
    /// @brief 
    /// @tparam scalar 
    /// @param direc from edge to point outside, normalized
    /// @param e0 
    /// @param e1 
    /// @param f0 face [f0, e0, e1]
    /// @param f1 face [f1, e1, e0]
    /// @param alpha 
    /// @return 
    template <typename scalar>
    inline scalar smooth_edge3_term(
        const Eigen::Ref<const Vector3<scalar>>& direc,
        const Eigen::Ref<const Vector3<scalar>>& e0,
        const Eigen::Ref<const Vector3<scalar>>& e1,
        const Eigen::Ref<const Vector3<scalar>>& f0,
        const Eigen::Ref<const Vector3<scalar>>& f1,
        const double alpha,
        const double beta)
    {
        const Vector3<scalar> t0 = point_line_closest_point_direction<scalar>(f0, e0, e1);
        const Vector3<scalar> t1 = point_line_closest_point_direction<scalar>(f1, e0, e1);
        scalar tangent_term = smooth_heaviside<scalar>(-direc.dot(t0) / t0.norm(), alpha, beta) *
                            smooth_heaviside<scalar>(-direc.dot(t1) / t1.norm(), alpha, beta);

        const Vector3<scalar> n0 = (e0 - f0).cross(e1 - f0);
        const Vector3<scalar> n1 = -(e0 - f1).cross(e1 - f1);
        scalar normal_term = smooth_heaviside<scalar>( (smooth_heaviside<scalar>(direc.dot(n0) / n0.norm(), alpha, beta) +
                             smooth_heaviside<scalar>(direc.dot(n1) / n1.norm(), alpha, beta) - 1), alpha, 0);

        return (e1 - e0).squaredNorm() * tangent_term * normal_term;
    }

    inline bool smooth_edge3_term_type(
        const Eigen::Ref<const Vector3<double>>& direc,
        const Eigen::Ref<const Vector3<double>>& e0,
        const Eigen::Ref<const Vector3<double>>& e1,
        const Eigen::Ref<const Vector3<double>>& f0,
        const Eigen::Ref<const Vector3<double>>& f1,
        const double alpha,
        const double beta)
    {
        const Vector3<double> t0 = point_line_closest_point_direction<double>(f0, e0, e1);
        const Vector3<double> t1 = point_line_closest_point_direction<double>(f1, e0, e1);
        if (-direc.dot(t0) / t0.norm() <= -alpha || -direc.dot(t1) / t1.norm() <= -alpha)
            return false;

        const Vector3<double> n0 = (e0 - f0).cross(e1 - f0);
        const Vector3<double> n1 = -(e0 - f1).cross(e1 - f1);
        if (smooth_heaviside<double>(direc.dot(n0) / n0.norm(), alpha, beta) + 
            smooth_heaviside<double>(direc.dot(n1) / n1.norm(), alpha, beta) <= 1 - alpha)
            return false;

        return true;
    }
}