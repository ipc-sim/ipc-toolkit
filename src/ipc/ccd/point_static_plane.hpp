#pragma once

#include <ipc/ccd/ccd.hpp>

namespace ipc {

bool point_static_plane_ccd(
    const VectorMax3d& p_t0,
    const VectorMax3d& p_t1,
    const VectorMax3d& plane_origin,
    const VectorMax3d& plane_normal,
    double& toi,
    const double conservative_rescaling = DEFAULT_CCD_CONSERVATIVE_RESCALING);

} // namespace ipc
