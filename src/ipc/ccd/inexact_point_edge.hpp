//
// NOTE: This method is provided for reference comparison and is not utilized by
// the high-level functionality. In compairson to Tight Inclusion CCD, this CCD
// method is not provably conservative and so can potentially produce false
// negatives (i.e., miss collisions) due to floating-point rounding error.
//

#pragma once

#include <Eigen/Core>

namespace ipc {

/// @brief Inexact continuous collision detection between a point and an edge in 2D.
/// @param[in] p_t0 The initial position of the point.
/// @param[in] e0_t0 The initial position of the first endpoint of the edge.
/// @param[in] e1_t0 The initial position of the second endpoint of the edge.
/// @param[in] p_t1 The final position of the point.
/// @param[in] e0_t1 The final position of the first endpoint of the edge.
/// @param[in] e1_t1 The final position of the second endpoint of the edge.
/// @param[out] toi Output time of impact.
/// @param[in] conservative_rescaling The conservative rescaling of the time of impact.
/// @return True if a collision was detected, false otherwise.
bool inexact_point_edge_ccd_2D(
    const Eigen::Vector2d& p_t0,
    const Eigen::Vector2d& e0_t0,
    const Eigen::Vector2d& e1_t0,
    const Eigen::Vector2d& p_t1,
    const Eigen::Vector2d& e0_t1,
    const Eigen::Vector2d& e1_t1,
    double& toi,
    const double conservative_rescaling);

} // namespace ipc
