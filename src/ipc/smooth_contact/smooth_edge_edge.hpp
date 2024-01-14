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
    /// assuming that edge [ea0, ea1] follows the orientation in face fa0 (so opposite orientation in fa1)
    /// edge [eb0, eb1] follows the orientation in face fb0 (so opposite orientation in fb1)
    /// @param fa0_normal, fa1_normal are the normals of two faces adjacent to edge [ea0, ea1]
    template <typename scalar>
    scalar smooth_edge_edge_potential_single_point(
        const Eigen::Ref<const Vector3<scalar>>& ea0,
        const Eigen::Ref<const Vector3<scalar>>& ea1,
        const Eigen::Ref<const Vector3<scalar>>& eb0,
        const Eigen::Ref<const Vector3<scalar>>& eb1,
        const Eigen::Ref<const Vector3<scalar>>& fa0_normal,
        const Eigen::Ref<const Vector3<scalar>>& fa1_normal,
        const Eigen::Ref<const Vector3<scalar>>& fb0_normal,
        const Eigen::Ref<const Vector3<scalar>>& fb1_normal,
        const ParameterType &params,
        EdgeEdgeDistanceType dtype)
    {
        constexpr double threshold_eps = 1e-3;

        const Vector3<scalar> u = ea1 - ea0;
        const Vector3<scalar> v = eb1 - eb0;
        const scalar a = u.squaredNorm();
        const scalar b = v.squaredNorm();
        const scalar mollifier_threshold = threshold_eps * a * b;
        const scalar cross_sqr_norm = u.cross(v).squaredNorm();
        const scalar mollifier_val = mollifier<scalar>(cross_sqr_norm / mollifier_threshold);
        Vector3<scalar> direc = edge_edge_closest_point_direction(ea0, ea1, eb0, eb1, dtype); // from edge a to edge b
        const scalar dist_sqr = edge_edge_sqr_distance(ea0, ea1, eb0, eb1, dtype); // get rid of me after verified!
        
        if constexpr (std::is_same<scalar, double>::value)
        {
            assert ((direc.squaredNorm() - dist_sqr) < 1e-12 * std::max(1., dist_sqr));
        }

        if (dist_sqr > params.eps)
            return scalar(0.);
        
        // vanishes if Phi < -alpha
        direc = direc / sqrt(dist_sqr);
        const scalar Phia0 = -direc.dot(fa0_normal.cross(u) / sqrt(a)) / params.alpha;
        const scalar Phia1 = -direc.dot(fa1_normal.cross(-u) / sqrt(a)) / params.alpha;
        const scalar Phib0 = direc.dot(fb0_normal.cross(v) / sqrt(b)) / params.alpha;
        const scalar Phib1 = direc.dot(fb1_normal.cross(-v) / sqrt(b)) / params.alpha;

        const scalar mollifier_a = mollifier<scalar>((point_edge_sqr_distance<scalar>(ea0, eb0, eb1) - dist_sqr) / a / threshold_eps) * mollifier<scalar>((point_edge_sqr_distance<scalar>(ea1, eb0, eb1) - dist_sqr) / a / threshold_eps);
        const scalar mollifier_b = mollifier<scalar>((point_edge_sqr_distance<scalar>(eb0, ea0, ea1) - dist_sqr) / b / threshold_eps) * mollifier<scalar>((point_edge_sqr_distance<scalar>(eb1, ea0, ea1) - dist_sqr) / b / threshold_eps);

        return 0.5 * u.norm() * v.norm() * inv_barrier<scalar>(dist_sqr / params.eps, params.r) * mollifier_val * (smooth_heaviside<scalar>(Phia0) * smooth_heaviside<scalar>(Phia1) * mollifier_a + smooth_heaviside<scalar>(Phib0) * smooth_heaviside<scalar>(Phib1) * mollifier_b);
    }
}