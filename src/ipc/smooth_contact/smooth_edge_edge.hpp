#pragma once

#include "edge.hpp"
#include "common.hpp"

namespace ipc {
    inline bool smooth_edge_edge_potential_type(
        const Eigen::Ref<const Vector3<double>>& ea0,
        const Eigen::Ref<const Vector3<double>>& ea1,
        const Eigen::Ref<const Vector3<double>>& eb0,
        const Eigen::Ref<const Vector3<double>>& eb1,
        const Eigen::Ref<const Vector3<double>>& fa0,
        const Eigen::Ref<const Vector3<double>>& fa1,
        const Eigen::Ref<const Vector3<double>>& fb0,
        const Eigen::Ref<const Vector3<double>>& fb1,
        const ParameterType &params,
        const EdgeEdgeDistanceType &dtype)
    {
        const double dist = sqrt(edge_edge_sqr_distance(ea0, ea1, eb0, eb1, dtype));
        if (dist >= params.dhat)
            return false;
        
        const Vector3<double> direc = edge_edge_closest_point_direction(ea0, ea1, eb0, eb1, dtype) / dist; // from edge a to edge b
        return smooth_edge3_term_type(direc, ea0, ea1, fa0, fa1, params.alpha, params.beta) && smooth_edge3_term_type(-direc, eb0, eb1, fb0, fb1, params.alpha, params.beta);
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
        const EdgeEdgeDistanceType &dtype)
    {
        const scalar dist_sqr = edge_edge_sqr_distance(ea0, ea1, eb0, eb1, dtype);
        const scalar barrier = inv_barrier<scalar>(sqrt(dist_sqr) / params.dhat, params.r);
        
        const Vector3<scalar> direc = edge_edge_closest_point_direction(ea0, ea1, eb0, eb1, dtype) / sqrt(dist_sqr); // from edge a to edge b
        const scalar out = smooth_edge3_term<scalar>(direc, ea0, ea1, fa0, fa1, params.alpha, params.beta) * 
                          smooth_edge3_term<scalar>(-direc, eb0, eb1, fb0, fb1, params.alpha, params.beta);

        const scalar mollifier_val = edge_edge_mollifier<scalar>(ea0, ea1, eb0, eb1, dist_sqr);

        // if constexpr (std::is_same<double,scalar>::value)
        // {
        //     if (dist_sqr < 1e-20)
        //     {
        //         logger().error("barrier {}, tangent {} {}, mollifier {}", barrier, smooth_edge3_term<scalar>(direc, ea0, ea1, fa0, fa1, params.alpha, params.beta, debug), smooth_edge3_term<scalar>(-direc, eb0, eb1, fb0, fb1, params.alpha, params.beta, debug), mollifier_val);
        //         logger().error("tangent types {} {}", smooth_edge3_term_type(direc, ea0, ea1, fa0, fa1, params.alpha, params.beta), smooth_edge3_term_type(-direc, eb0, eb1, fb0, fb1, params.alpha, params.beta));
        //         std::cout << ea0.transpose() << "\n"  << ea1.transpose() << "\n"  << eb0.transpose() << "\n"  << eb1.transpose() << "\n"  << fa0.transpose() << "\n"  << fa1.transpose() << "\n"  << fb0.transpose() << "\n"  << fb1.transpose() << "\n";
        //         std::cout << direc.transpose() << "\n" <<
        //                      point_line_closest_point_direction<double>(fa0, ea0, ea1).transpose() << "\n" <<
        //                     point_line_closest_point_direction<double>(fa1, ea0, ea1).transpose() << "\n" <<
        //                     point_line_closest_point_direction<double>(fb0, eb0, eb1).transpose() << "\n" <<
        //                     point_line_closest_point_direction<double>(fb1, eb0, eb1).transpose() << "\n";
        //     }
        // }
        
        return barrier * out * mollifier_val;
    }
}