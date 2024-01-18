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

    template <typename scalar>
    scalar smooth_edge_edge_potential_tangent_term(
        const Eigen::Ref<const Vector3<scalar>>& f0,
        const Eigen::Ref<const Vector3<scalar>>& f1,
        const Eigen::Ref<const Vector3<scalar>>& f2,
        const Eigen::Ref<const Vector3<scalar>>& direc,
        const double &alpha,
        const HEAVISIDE_TYPE &type)
    {
        switch (type)
        {
        case HEAVISIDE_TYPE::ZERO:
            return scalar(0.);
        case HEAVISIDE_TYPE::ONE:
            return scalar(1.);
        case HEAVISIDE_TYPE::VARIANT:
        default:
            Vector3<scalar> t = point_line_closest_point_direction<scalar>(f0, f1, f2);
            return smooth_heaviside<scalar>(direc.dot(t) / t.norm() / alpha);
        }
    }

    template <typename scalar>
    scalar smooth_edge_edge_potential_normal_term(
        const Eigen::Ref<const Vector3<scalar>>& f0,
        const Eigen::Ref<const Vector3<scalar>>& f1,
        const Eigen::Ref<const Vector3<scalar>>& f2,
        const Eigen::Ref<const Vector3<scalar>>& direc,
        const double &alpha,
        const HEAVISIDE_TYPE &type)
    {
        switch (type)
        {
        case HEAVISIDE_TYPE::ZERO:
            return scalar(0.);
        case HEAVISIDE_TYPE::ONE:
            return scalar(1.);
        case HEAVISIDE_TYPE::VARIANT:
        default:
            const Vector3<scalar> normal = (f1 - f0).cross(f2 - f0);
            return smooth_heaviside<scalar>(direc.dot(normal) / normal.norm() / alpha);
        }
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
        const std::array<double, 2> &dhats,
        const EdgeEdgeDistanceType &dtype,
        const std::array<PointEdgeDistanceType, 4> &edge_dtypes,
        const std::array<HEAVISIDE_TYPE, 4> &tangent_types,
        const std::array<HEAVISIDE_TYPE, 4> &normal_types)
    {
        const Vector3<scalar> u = ea1 - ea0;
        const Vector3<scalar> v = eb1 - eb0;
        const scalar a = u.squaredNorm();
        const scalar b = v.squaredNorm();
        const scalar dist_sqr = edge_edge_sqr_distance(ea0, ea1, eb0, eb1, dtype);
        
        // vanishes if Phi < -alpha
        const Vector3<scalar> direc = edge_edge_closest_point_direction(ea0, ea1, eb0, eb1, dtype) / sqrt(dist_sqr); // from edge a to edge b
        scalar out = scalar(0.);
        if (tangent_types[0] != HEAVISIDE_TYPE::ZERO && tangent_types[1] != HEAVISIDE_TYPE::ZERO)
            out += smooth_edge_edge_potential_tangent_term<scalar>(fa0, ea0, ea1, -direc, params.alpha, tangent_types[0]) *
                    smooth_edge_edge_potential_tangent_term<scalar>(fa1, ea1, ea0, -direc, params.alpha, tangent_types[1]) *
                    inv_barrier<scalar>(dist_sqr / intpow(dhats[0], 2), params.r);
        if (tangent_types[2] != HEAVISIDE_TYPE::ZERO && tangent_types[3] != HEAVISIDE_TYPE::ZERO)
            out += smooth_edge_edge_potential_tangent_term<scalar>(fb0, eb0, eb1, direc, params.alpha, tangent_types[2]) *
                    smooth_edge_edge_potential_tangent_term<scalar>(fb1, eb1, eb0, direc, params.alpha, tangent_types[3]) *
                    inv_barrier<scalar>(dist_sqr / intpow(dhats[1], 2), params.r);

        const scalar normal_penalty = (smooth_edge_edge_potential_normal_term<scalar>(fa0, ea0, ea1, direc, params.alpha, normal_types[0]) + 
                                    smooth_edge_edge_potential_normal_term<scalar>(fa1, ea1, ea0, direc, params.alpha, normal_types[1])) * 
                                    (smooth_edge_edge_potential_normal_term<scalar>(fb0, eb0, eb1, -direc, params.alpha, normal_types[2]) + 
                                    smooth_edge_edge_potential_normal_term<scalar>(fb1, eb1, eb0, -direc, params.alpha, normal_types[3]));

        const scalar mollifier_a = mollifier<scalar>((point_edge_sqr_distance<scalar>(ea0, eb0, eb1, edge_dtypes[0]) - dist_sqr) / a / mollifier_threshold_eps) * mollifier<scalar>((point_edge_sqr_distance<scalar>(ea1, eb0, eb1, edge_dtypes[1]) - dist_sqr) / a / mollifier_threshold_eps);
        const scalar mollifier_b = mollifier<scalar>((point_edge_sqr_distance<scalar>(eb0, ea0, ea1, edge_dtypes[2]) - dist_sqr) / b / mollifier_threshold_eps) * mollifier<scalar>((point_edge_sqr_distance<scalar>(eb1, ea0, ea1, edge_dtypes[3]) - dist_sqr) / b / mollifier_threshold_eps);

        return 0.5 * sqrt(a * b) * out * normal_penalty * mollifier_a * mollifier_b;
    }
}