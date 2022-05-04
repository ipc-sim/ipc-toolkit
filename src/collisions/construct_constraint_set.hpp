#pragma once

#include <Eigen/Core>

#include <tbb/enumerable_thread_specific.h>

#include <ipc/collision_mesh.hpp>
#include <ipc/collisions/constraints.hpp>
#include <ipc/broad_phase/broad_phase.hpp>

namespace ipc {

/// @brief Construct a set of constraints used to compute the barrier potential.
/// @param[in] mesh The collision mesh.
/// @param[in] V Vertices of the collision mesh.
/// @param[in] dhat The activation distance of the barrier.
/// @param[out] constraint_set The constructed set of constraints (any existing constraints will be cleared).
/// @param[in] dmin Minimum distance.
/// @param[in] method Broad-phase method to use.
void construct_constraint_set(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const double dhat,
    Constraints& constraint_set,
    const double dmin = 0,
    const BroadPhaseMethod method = BroadPhaseMethod::HASH_GRID);

/// @brief Construct a set of constraints used to compute the barrier potential.
/// @param[in] candidates Distance candidates from which the constraint set is built.
/// @param[in] mesh The collision mesh.
/// @param[in] V Vertices of the collision mesh.
/// @param[in] dhat The activation distance of the barrier.
/// @param[out] constraint_set The constructed set of constraints (any existing constraints will be cleared).
/// @param[in]  dmin  Minimum distance.
void construct_constraint_set(
    const Candidates& candidates,
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const double dhat,
    Constraints& constraint_set,
    const double dmin = 0);

///////////////////////////////////////////////////////////////////////////////

struct ConstraintBuilder {
    // Store the indices to VV and EV pairs to avoid duplicates.
    unordered_map<VertexVertexConstraint, long> vv_to_index;
    unordered_map<EdgeVertexConstraint, long> ev_to_index;
    Constraints constraint_set;
};

void edge_vertex_candiates_to_constraints(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const std::vector<EdgeVertexCandidate>& candidates,
    const std::function<bool(double)>& is_active,
    const size_t start_i,
    const size_t end_i,
    ConstraintBuilder& constraint_builder);

void edge_edge_candiates_to_constraints(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const std::vector<EdgeEdgeCandidate>& candidates,
    const std::function<bool(double)>& is_active,
    const size_t start_i,
    const size_t end_i,
    ConstraintBuilder& constraint_builder);

void face_vertex_candiates_to_constraints(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const std::vector<FaceVertexCandidate>& candidates,
    const std::function<bool(double)>& is_active,
    const size_t start_i,
    const size_t end_i,
    ConstraintBuilder& constraint_builder);

void merge_thread_local_constraints(
    const tbb::enumerable_thread_specific<ConstraintBuilder>& local_storage,
    Constraints& constraints);

} // namespace ipc
