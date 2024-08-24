#pragma once

#include <ipc/utils/logger.hpp>

namespace ipc {

inline bool check_initial_distance(
    const double initial_distance, const double min_distance, double& toi)
{
    if (initial_distance > min_distance) {
        return false;
    }

    logger().warn(
        "Initial distance {} ≤ d_min={}, returning toi=0!", initial_distance,
        min_distance);

    toi = 0; // Initially touching

    return true;
}

} // namespace ipc