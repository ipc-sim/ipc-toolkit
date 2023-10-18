//
// NOTE: This method is provided for reference comparison and is not utilized by
// the high-level functionality. In compairson to Tight Inclusion CCD, this CCD
// method is not provably conservative and so can potentially produce false
// negatives (i.e., miss collisions) due to floating-point rounding error.
//

#pragma once

#include <Eigen/Core>

namespace ipc {

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
