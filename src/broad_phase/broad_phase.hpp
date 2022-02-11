#pragma once

#include <Eigen/Core>

#include <ipc/broad_phase/collision_candidate.hpp>
#include <ipc/collision_mesh.hpp>

namespace ipc {

enum class BroadPhaseMethod {
    BRUTE_FORCE,
    HASH_GRID,
    SPATIAL_HASH,
#ifdef IPC_TOOLKIT_WITH_CUDA
    SWEEP_AND_TINIEST_QUEUE,
#endif
};

/// @brief Construct a set of discrete collision detection candidates.
/// @param[in] mesh The surface of the contact mesh.
/// @param[in] V Surface Vertex positions at start as rows of a matrix.
/// @param[out] canidates The constructed candidate set as output.
/// @param[in] inflation_radius Amount to inflate the bounding boxes.
/// @param[in] method Broad phase method to use.
void construct_collision_candidates(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    Candidates& candidates,
    double inflation_radius = 0,
    const BroadPhaseMethod& method = BroadPhaseMethod::HASH_GRID);

/// @brief Construct a set of continous collision detection candidates.
/// @note Assumes the trajectory is linear.
/// @param[in] mesh The surface of the contact mesh.
/// @param[in] V0 Surface vertex positions at start as rows of a matrix.
/// @param[in] V1 Surface vertex positions at end as rows of a matrix.
/// @param[out] canidates The constructed candidate set as output.
/// @param[in] inflation_radius Amount to inflate the bounding boxes.
/// @param[in] method Broad phase method to use.
void construct_collision_candidates(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    Candidates& candidates,
    double inflation_radius = 0,
    const BroadPhaseMethod& method = BroadPhaseMethod::HASH_GRID);

} // namespace ipc
