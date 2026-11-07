#pragma once

#include <ipc/utils/logger.hpp>

#include <cmath>

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

/// @brief Shrink a minimum separation until a positive initial distance clears it.
///
/// A component query that cannot report toi = 0 trades separation margin for
/// a positive time of impact: we halve the minimum separation until the
/// initial distance exceeds it (mirroring Tight Inclusion's small-TOI rerun).
///
/// @param initial_distance The initial (unsquared) distance between the objects.
/// @param min_distance The requested minimum separation.
/// @return The minimum separation, halved just enough to fall below a positive initial_distance (unchanged otherwise).
inline double
shrink_min_distance_below(const double initial_distance, double min_distance)
{
    if (initial_distance <= 0 || initial_distance > min_distance) {
        return min_distance;
    }
    // min_distance / 2^k < initial_distance for the smallest such k ≥ 1.
    const int k =
        static_cast<int>(std::floor(std::log2(min_distance / initial_distance)))
        + 1;
    return std::ldexp(min_distance, -k); // halve k times
}

} // namespace ipc