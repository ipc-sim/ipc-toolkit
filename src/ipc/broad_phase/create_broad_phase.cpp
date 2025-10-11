#include "create_broad_phase.hpp"

#include <ipc/broad_phase/brute_force.hpp>
#include <ipc/broad_phase/bvh.hpp>
#include <ipc/broad_phase/hash_grid.hpp>
#include <ipc/broad_phase/spatial_hash.hpp>
#include <ipc/broad_phase/sweep_and_prune.hpp>
#include <ipc/broad_phase/sweep_and_tiniest_queue.hpp>

namespace ipc {

std::shared_ptr<BroadPhase>
create_broad_phase(const BroadPhaseMethod& broad_phase_method)
{
    switch (broad_phase_method) {
    case BroadPhaseMethod::BRUTE_FORCE:
        return std::make_shared<BruteForce>();
    case BroadPhaseMethod::HASH_GRID:
        return std::make_shared<HashGrid>();
    case BroadPhaseMethod::SPATIAL_HASH:
        return std::make_shared<SpatialHash>();
    case BroadPhaseMethod::BVH:
        return std::make_shared<BVH>();
    case BroadPhaseMethod::SWEEP_AND_PRUNE:
        return std::make_shared<SweepAndPrune>();
#ifdef IPC_TOOLKIT_WITH_CUDA
    case BroadPhaseMethod::SWEEP_AND_TINIEST_QUEUE:
        return std::make_shared<SweepAndTiniestQueue>();
#endif
    default:
        log_and_throw_error("Unknown broad phase type!");
    }
}

} // namespace ipc
