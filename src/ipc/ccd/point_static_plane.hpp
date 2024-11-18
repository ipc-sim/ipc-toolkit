#pragma once

#include <ipc/ccd/tight_inclusion_ccd.hpp>

namespace ipc {

/// @brief Computes the time of impact between a point and a static plane in 3D using continuous collision detection.
/// @param[in] p_t0 The initial position of the point.
/// @param[in] p_t1 The final position of the point.
/// @param[in] plane_origin The origin of the plane.
/// @param[in] plane_normal The normal of the plane.
/// @param[out] toi Output time of impact.
/// @param[in] conservative_rescaling Conservative rescaling of the time of impact.
/// @return True if a collision was detected, false otherwise.
bool point_static_plane_ccd(
    const VectorMax3d& p_t0,
    const VectorMax3d& p_t1,
    const VectorMax3d& plane_origin,
    const VectorMax3d& plane_normal,
    double& toi,
    const double conservative_rescaling =
        TightInclusionCCD::DEFAULT_CONSERVATIVE_RESCALING);

} // namespace ipc
