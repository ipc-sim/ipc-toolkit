#pragma once

#include <ipc/broad_phase/broad_phase.hpp>
#include <ipc/collision_mesh.hpp>

#include <Eigen/Core>
#include <Eigen/Sparse>

/// Incremental Potential Contact functions
namespace ipc {

// ============================================================================
// Collision detection

/// @brief Determine if the step is collision free.
/// @note Assumes the trajectory is linear.
/// @param mesh The collision mesh.
/// @param vertices_t0 Surface vertex vertices at start as rows of a matrix.
/// @param vertices_t1 Surface vertex vertices at end as rows of a matrix.
/// @param broad_phase_method The broad phase method to use.
/// @param min_distance The minimum distance allowable between any two elements.
/// @param tolerance The tolerance for the CCD algorithm.
/// @param max_iterations The maximum number of iterations for the CCD algorithm.
/// @returns True if <b>any</b> collisions occur.
bool is_step_collision_free(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const BroadPhaseMethod broad_phase_method = DEFAULT_BROAD_PHASE_METHOD,
    const double min_distance = 0.0,
    const double tolerance = DEFAULT_CCD_TOLERANCE,
    const long max_iterations = DEFAULT_CCD_MAX_ITERATIONS);

/// @brief Computes a maximal step size that is collision free.
/// @note Assumes the trajectory is linear.
/// @param mesh The collision mesh.
/// @param vertices_t0 Vertex vertices at start as rows of a matrix. Assumes vertices_t0 is intersection free.
/// @param vertices_t1 Surface vertex vertices at end as rows of a matrix.
/// @param broad_phase_method The broad phase method to use.
/// @param min_distance The minimum distance allowable between any two elements.
/// @param tolerance The tolerance for the CCD algorithm.
/// @param max_iterations The maximum number of iterations for the CCD algorithm.
/// @returns A step-size \f$\in [0, 1]\f$ that is collision free. A value of 1.0 if a full step and 0.0 is no step.
double compute_collision_free_stepsize(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const BroadPhaseMethod broad_phase_method = DEFAULT_BROAD_PHASE_METHOD,
    const double min_distance = 0.0,
    const double tolerance = DEFAULT_CCD_TOLERANCE,
    const long max_iterations = DEFAULT_CCD_MAX_ITERATIONS);

// ============================================================================
// Utilities

/// @brief Determine if the mesh has self intersections.
/// @param mesh The collision mesh.
/// @param vertices Vertices of the collision mesh.
/// @param broad_phase_method The broad phase method to use.
/// @return A boolean for if the mesh has intersections.
bool has_intersections(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const BroadPhaseMethod broad_phase_method = DEFAULT_BROAD_PHASE_METHOD);

} // namespace ipc
