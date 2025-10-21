#pragma once

#include <ipc/config.hpp>
#include <ipc/broad_phase/broad_phase.hpp>

namespace ipc {

enum class BroadPhaseMethod : uint8_t {
    BRUTE_FORCE,
    HASH_GRID,
    SPATIAL_HASH,
    BVH,
    SWEEP_AND_PRUNE,
    SWEEP_AND_TINIEST_QUEUE
};

std::shared_ptr<BroadPhase>
create_broad_phase(const BroadPhaseMethod& broad_phase_method);

} // namespace ipc
