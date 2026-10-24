#pragma once

#include <ipc/config.hpp>
#include <ipc/broad_phase/broad_phase.hpp>

namespace ipc {

/// @brief The broad phase methods create_broad_phase() can construct.
enum class BroadPhaseMethod : uint8_t {
    BRUTE_FORCE,
    HASH_GRID,
    SPATIAL_HASH,
    LBVH,
    SWEEP_AND_PRUNE,
    SWEEP_AND_TINIEST_QUEUE, ///< Requires CUDA (IPC_TOOLKIT_WITH_CUDA).
    LBVH_CUDA,               ///< Requires CUDA (IPC_TOOLKIT_WITH_CUDA).
    /// @brief The number of methods; not a method itself.
    NUM_BROAD_PHASE_METHODS
};

/// @brief Construct a broad phase of the given method.
/// @param broad_phase_method The method to construct.
/// @return The broad phase.
/// @throws std::runtime_error if the method requires CUDA and the library was
/// built without it, or if the method is not a valid BroadPhaseMethod.
std::shared_ptr<BroadPhase>
create_broad_phase(const BroadPhaseMethod& broad_phase_method);

} // namespace ipc
