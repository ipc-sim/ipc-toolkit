#pragma once

#include "smooth_point_face.hpp"
#include <ipc/utils/distance_autodiff.hpp>

namespace ipc {
    /// @brief Old IPC formulation. Compute edge-edge potential for edge [ea0, ea1] and [eb0, eb1]
    /// @param fa0_normal, fa1_normal are the normals of two faces adjacent to edge [ea0, ea1]
    template <typename scalar>
    scalar edge_edge_potential_single_point(
        const Eigen::Ref<const Vector3<scalar>>& ea0,
        const Eigen::Ref<const Vector3<scalar>>& ea1,
        const Eigen::Ref<const Vector3<scalar>>& eb0,
        const Eigen::Ref<const Vector3<scalar>>& eb1,
        const ParameterType &params,
        EdgeEdgeDistanceType dtype)
    {
        const Vector3<scalar> u = ea1 - ea0;
        const Vector3<scalar> v = eb1 - eb0;
        const scalar mollifier_threshold = 1e-3 * u.squaredNorm() * v.squaredNorm();
        const scalar cross_sqr_norm = u.cross(v).squaredNorm();
        const scalar mollifier_val = mollifier<scalar>(cross_sqr_norm / mollifier_threshold);
        const scalar dist_sqr = edge_edge_sqr_distance(ea0, ea1, eb0, eb1, dtype);

        return inv_barrier<scalar>(dist_sqr / params.eps, params.r) * mollifier_val;
    }

    /// @brief Compute edge-edge potential for edge [ea0, ea1] and [eb0, eb1]
    /// edge [ea0, ea1] follows the orientation in face fa0 = [fa0, ea0, ea1], so opposite orientation in fa1 = [fa1, ea1, ea0]
    /// edge [eb0, eb1] follows the orientation in face fb0 = [fb0, eb0, eb1], so opposite orientation in fb1 = [fb1, eb1, eb0]
    template <typename scalar>
    scalar smooth_edge_edge_potential_single_point(
        const Eigen::Ref<const Vector3<scalar>>& ea0,
        const Eigen::Ref<const Vector3<scalar>>& ea1,
        const Eigen::Ref<const Vector3<scalar>>& eb0,
        const Eigen::Ref<const Vector3<scalar>>& eb1,
        const Eigen::Ref<const Vector3<scalar>>& fa0,
        const Eigen::Ref<const Vector3<scalar>>& fa1,
        const Eigen::Ref<const Vector3<scalar>>& fb0,
        const Eigen::Ref<const Vector3<scalar>>& fb1,
        const ParameterType &params,
        const EdgeEdgeDistanceType &dtype,
        const std::array<PointEdgeDistanceType, 4> &edge_dtypes)
    {
        constexpr double threshold_eps = 1e-2;

        const Vector3<scalar> u = ea1 - ea0;
        const Vector3<scalar> v = eb1 - eb0;
        const scalar a = u.squaredNorm();
        const scalar b = v.squaredNorm();
        const scalar dist_sqr = edge_edge_sqr_distance(ea0, ea1, eb0, eb1, dtype);

        if (dist_sqr > params.eps)
            return scalar(0.);
        
        // vanishes if Phi < -alpha
        const Vector3<scalar> direc = edge_edge_closest_point_direction(ea0, ea1, eb0, eb1, dtype) / sqrt(dist_sqr); // from edge a to edge b
        scalar tangent_penalty = scalar(1.);
        {
            const scalar Phia0 = -direc.dot(point_line_closest_point_direction<scalar>(fa0, ea0, ea1).normalized()) / params.alpha;
            const scalar Phia1 = -direc.dot(point_line_closest_point_direction<scalar>(fa1, ea0, ea1).normalized()) / params.alpha;
            const scalar Phib0 = direc.dot(point_line_closest_point_direction<scalar>(fb0, eb0, eb1).normalized()) / params.alpha;
            const scalar Phib1 = direc.dot(point_line_closest_point_direction<scalar>(fb1, eb0, eb1).normalized()) / params.alpha;
            tangent_penalty = smooth_heaviside<scalar>(Phia0) * smooth_heaviside<scalar>(Phia1) + smooth_heaviside<scalar>(Phib0) * smooth_heaviside<scalar>(Phib1);
        }

        scalar normal_penalty = scalar(1.);
        {
            const Vector3<scalar> n_a0 = (ea0 - fa0).cross(ea1 - fa0).normalized();
            const Vector3<scalar> n_a1 = (ea1 - fa1).cross(ea0 - fa1).normalized();
            const Vector3<scalar> n_b0 = (eb0 - fb0).cross(eb1 - fb0).normalized();
            const Vector3<scalar> n_b1 = (eb1 - fb1).cross(eb0 - fb1).normalized();
            normal_penalty = smooth_heaviside<scalar>(direc.dot(n_a0)) * smooth_heaviside<scalar>(direc.dot(n_a1)) * smooth_heaviside<scalar>(-direc.dot(n_b0)) * smooth_heaviside<scalar>(-direc.dot(n_b1));
        }

        const scalar mollifier_a = mollifier<scalar>((point_edge_sqr_distance<scalar>(ea0, eb0, eb1, edge_dtypes[0]) - dist_sqr) / a / threshold_eps) * mollifier<scalar>((point_edge_sqr_distance<scalar>(ea1, eb0, eb1, edge_dtypes[1]) - dist_sqr) / a / threshold_eps);
        const scalar mollifier_b = mollifier<scalar>((point_edge_sqr_distance<scalar>(eb0, ea0, ea1, edge_dtypes[2]) - dist_sqr) / b / threshold_eps) * mollifier<scalar>((point_edge_sqr_distance<scalar>(eb1, ea0, ea1, edge_dtypes[3]) - dist_sqr) / b / threshold_eps);

        return 0.5 * sqrt(a * b) * inv_barrier<scalar>(dist_sqr / params.eps, params.r) * tangent_penalty * normal_penalty * mollifier_a * mollifier_b;
    }
}