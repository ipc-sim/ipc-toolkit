#pragma once

#include <ipc/collision_mesh.hpp>
#include <ipc/broad_phase/default_broad_phase.hpp>
#include <ipc/ccd/default_narrow_phase_ccd.hpp>

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
/// @param min_distance The minimum distance allowable between any two elements.
/// @param broad_phase_method The broad phase method to use.
/// @param narrow_phase_ccd The narrow phase CCD algorithm to use.
/// @returns True if <b>any</b> collisions occur.
bool is_step_collision_free(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const double min_distance = 0.0,
    const std::shared_ptr<BroadPhase> broad_phase = make_default_broad_phase(),
    const NarrowPhaseCCD& narrow_phase_ccd = DEFAULT_NARROW_PHASE_CCD);

/// @brief Computes a maximal step size that is collision free.
/// @note Assumes the trajectory is linear.
/// @param mesh The collision mesh.
/// @param vertices_t0 Vertex vertices at start as rows of a matrix. Assumes vertices_t0 is intersection free.
/// @param vertices_t1 Surface vertex vertices at end as rows of a matrix.
/// @param min_distance The minimum distance allowable between any two elements.
/// @param broad_phase_method The broad phase method to use.
/// @param narrow_phase_ccd The narrow phase CCD algorithm to use.
/// @returns A step-size \f$\in [0, 1]\f$ that is collision free. A value of 1.0 if a full step and 0.0 is no step.
double compute_collision_free_stepsize(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const double min_distance = 0.0,
    const std::shared_ptr<BroadPhase> broad_phase = make_default_broad_phase(),
    const NarrowPhaseCCD& narrow_phase_ccd = DEFAULT_NARROW_PHASE_CCD);

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
    const std::shared_ptr<BroadPhase> broad_phase = make_default_broad_phase());

} // namespace ipc
