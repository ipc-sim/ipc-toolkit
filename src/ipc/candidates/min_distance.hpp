#pragma once

#include <ipc/collision_mesh.hpp>
#include <ipc/utils/logger.hpp>

#include <mutex>
#include <optional>

namespace ipc {

/// @brief Warn (at most once per process) that min_distance is deprecated.
inline void warn_min_distance_deprecated()
{
    static std::once_flag flag;
    std::call_once(flag, [] {
        logger().warn(
            "The min_distance parameter is deprecated: give the collision mesh "
            "a thickness with CollisionMesh::set_thickness() instead. A "
            "min_distance overrides the mesh's thickness.");
    });
}

/// @brief Warn (at most once per process) that dmin is deprecated.
inline void warn_dmin_deprecated()
{
    static std::once_flag flag;
    std::call_once(flag, [] {
        logger().warn(
            "The dmin parameter is deprecated: give the collision mesh a "
            "thickness with CollisionMesh::set_thickness() instead. A dmin "
            "overrides the mesh's thickness.");
    });
}

/// @brief The distance a pair of primitives must keep.
/// @param min_distance The deprecated uniform minimum distance, which overrides the mesh's thickness when given.
/// @param candidate The pair of primitives.
/// @param mesh The collision mesh whose primitives the candidate indexes.
/// @return min_distance if given, else the pair's thickness (zero without one).
template <typename Candidate>
double pair_min_distance(
    const std::optional<double>& min_distance,
    const Candidate& candidate,
    const CollisionMesh& mesh)
{
    return min_distance.has_value() ? *min_distance
                                    : mesh.pair_thickness(candidate);
}

} // namespace ipc
