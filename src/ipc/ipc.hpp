#pragma once

// NOTE: Include this so the user can just include ipc.hpp
#include <ipc/collisions/constraints.hpp>
#include <ipc/friction/constraints.hpp>

#include <ipc/broad_phase/broad_phase.hpp>
#include <ipc/collision_mesh.hpp>

#include <Eigen/Core>
#include <Eigen/Sparse>

/// Incremental Potential Contact functions
namespace ipc {

///////////////////////////////////////////////////////////////////////////////
// Collision detection

/// @brief Determine if the step is collision free.
/// @note Assumes the trajectory is linear.
/// @param[in] mesh The collision mesh.
/// @param[in] V0 Surface vertex positions at start as rows of a matrix.
/// @param[in] V1 Surface vertex positions at end as rows of a matrix.
/// @returns True if <b>any</b> collisions occur.
bool is_step_collision_free(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const BroadPhaseMethod method = BroadPhaseMethod::HASH_GRID,
    const double tolerance = 1e-6,
    const long max_iterations = 1e7);

/// @brief Determine if the step is collision free from a set of candidates.
/// @note Assumes the trajectory is linear.
/// @param[in] candidates Set of candidates to check for collisions.
/// @param[in] mesh The collision mesh.
/// @param[in] V0 Surface vertex positions at start as rows of a matrix.
/// @param[in] V1 Surface vertex positions at end as rows of a matrix.
/// @returns True if <b>any</b> collisions occur.
bool is_step_collision_free(
    const Candidates& candidates,
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const double tolerance = 1e-6,
    const long max_iterations = 1e7);

/// @brief Computes a maximal step size that is collision free.
/// @note Assumes the trajectory is linear.
/// @param[in] mesh The collision mesh.
/// @param[in] V0 Vertex positions at start as rows of a matrix. Assumes V0 is intersection free.
/// @param[in] V1 Surface vertex positions at end as rows of a matrix.
/// @returns A step-size \f$\in [0, 1]\f$ that is collision free. A value of 1.0 if a full step and 0.0 is no step.
double compute_collision_free_stepsize(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const BroadPhaseMethod method = BroadPhaseMethod::HASH_GRID,
    const double tolerance = 1e-6,
    const long max_iterations = 1e7);

/// @brief Computes a maximal step size that is collision free using a set of collision candidates.
/// @note Assumes the trajectory is linear.
/// @param[in] candidates Set of candidates to check for collisions.
/// @param[in] mesh The collision mesh.
/// @param[in] V0 Vertex positions at start as rows of a matrix. Assumes V0 is intersection free.
/// @param[in] V1 Surface vertex positions at end as rows of a matrix.
/// @returns A step-size \f$\in [0, 1]\f$ that is collision free. A value of 1.0 if a full step and 0.0 is no step.
double compute_collision_free_stepsize(
    const Candidates& candidates,
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const double tolerance = 1e-6,
    const long max_iterations = 1e7);

///////////////////////////////////////////////////////////////////////////////
// Utilities

/// @brief Determine if the mesh has self intersections.
/// @param[in] mesh The collision mesh.
/// @param[in] V Vertices of the collision mesh.
/// @return A boolean for if the mesh has intersections.
bool has_intersections(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const BroadPhaseMethod method = BroadPhaseMethod::HASH_GRID);

} // namespace ipc
