#include "create_broad_phase.hpp"

#include <ipc/broad_phase/brute_force.hpp>
#include <ipc/broad_phase/bvh.hpp>
#include <ipc/broad_phase/hash_grid.hpp>
#include <ipc/broad_phase/lbvh.hpp>
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
    case BroadPhaseMethod::LBVH:
        return std::make_shared<LBVH>();
    case BroadPhaseMethod::SWEEP_AND_PRUNE:
        return std::make_shared<SweepAndPrune>();
    case BroadPhaseMethod::SWEEP_AND_TINIEST_QUEUE:
#ifdef IPC_TOOLKIT_WITH_CUDA
        return std::make_shared<SweepAndTiniestQueue>();
#else
        log_and_throw_error(
            "Sweep and Tiniest Queue broad phase requires CUDA! Enable it through CMake option IPC_TOOLKIT_WITH_CUDA.");
#endif
    default:
        log_and_throw_error("Unknown broad phase type!");
    }
}

} // namespace ipc
