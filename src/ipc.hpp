#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <ipc/collision_constraint.hpp>
// NOTE: Include this so the user can just include ipc.hpp
#include <ipc/friction/friction.hpp>
#include <ipc/broad_phase/broad_phase.hpp>
#include <ipc/collision_mesh.hpp>

/// Incremental Potential Contact functions
namespace ipc {

/// @brief Construct a set of constraints used to compute the barrier potential.
/// @param[in] mesh The surface of the contact mesh.
/// @param[in] V Surface vertex positions as rows of a matrix.
/// @param[in] dhat The activation distance of the barrier.
/// @param[out] constraint_set
///     The constructed set of constraints (any existing constraints will be
///     cleared).
/// @param[in] dmin Minimum distance.
/// @param[in] method Broad-phase method to use.
void construct_constraint_set(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const double dhat,
    Constraints& constraint_set,
    const double dmin = 0,
    const BroadPhaseMethod& method = BroadPhaseMethod::HASH_GRID);

/// @brief Construct a set of constraints used to compute the barrier potential.
/// @param[in] candidates Distance candidates from which the constraint set is
///                       built.
/// @param[in] mesh The surface of the contact mesh.
/// @param[in] V Surface vertex positions as rows of a matrix.
/// @param[in] dhat The activation distance of the barrier.
/// @param[out] constraint_set
///     The constructed set of constraints (any existing constraints will be
///     cleared).
/// @param[in]  dmin  Minimum distance.
void construct_constraint_set(
    const Candidates& candidates,
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const double dhat,
    Constraints& constraint_set,
    const double dmin = 0);

/// @brief Compute the barrier potential for a given constraint set.
///
/// @param[in] mesh The surface of the contact mesh.
/// @param[in] V Surface vertex positions as rows of a matrix.
/// @param[in] constraint_set The set of constraints.
/// @param[in] dhat The activation distance of the barrier.
/// @returns The sum of all barrier potentials (not scaled by the barrier
/// stiffness).
double compute_barrier_potential(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const Constraints& constraint_set,
    const double dhat);

/// @brief Compute the gradient of the barrier potential.
///
/// @param[in] mesh The surface of the contact mesh.
/// @param[in] V Surface vertex positions as rows of a matrix.
/// @param[in] constraint_set The set of constraints.
/// @param[in] dhat The activation distance of the barrier.
/// @returns The gradient of all barrier potentials (not scaled by the barrier
/// stiffness). This will have a size of |V|.
Eigen::VectorXd compute_barrier_potential_gradient(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const Constraints& constraint_set,
    const double dhat);

/// @brief Compute the hessian of the barrier potential.
///
/// @param[in] mesh The surface of the contact mesh.
/// @param[in] V Surface vertex positions as rows of a matrix.
/// @param[in] constraint_set The set of constraints.
/// @param[in] dhat The activation distance of the barrier.
/// @param[in] project_hessian_to_psd Make sure the hessian is positive
///            semi-definite.
/// @returns The hessian of all barrier potentials (not scaled by the barrier
///     stiffness). This will have a size of |V|x|V|.
Eigen::SparseMatrix<double> compute_barrier_potential_hessian(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const Constraints& constraint_set,
    const double dhat,
    const bool project_hessian_to_psd = true);

///////////////////////////////////////////////////////////////////////////////
// Collision detection

/// @brief Determine if the step is collision free.
///
/// V* can either be the surface vertices or the entire mesh vertices.
///
/// @note Assumes the trajectory is linear.
///
/// @param[in] mesh The surface of the contact mesh.
/// @param[in] V0 Surface vertex positions at start as rows of a matrix.
/// @param[in] V1 Surface vertex positions at end as rows of a matrix.
/// @returns True if <b>any</b> collisions occur.
bool is_step_collision_free(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const BroadPhaseMethod& method = BroadPhaseMethod::HASH_GRID,
    const double tolerance = 1e-6,
    const long max_iterations = 1e7);

/// @brief Determine if the step is collision free from a set of candidates.
///
/// V can either be the surface vertices or the entire mesh vertices. The edges
/// and face should be only for the surface elements.
///
/// @note Assumes the trajectory is linear.
///
/// @param[in] candidates Set of candidates to check for collisions.
/// @param[in] mesh The surface of the contact mesh.
/// @param[in] V0 Surface vertex positions at start as rows of a matrix.
/// @param[in] V1 Surface vertex positions at end as rows of a matrix.
///
/// @returns True if <b>any</b> collisions occur.
bool is_step_collision_free(
    const Candidates& candidates,
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const double tolerance = 1e-6,
    const long max_iterations = 1e7);

/// @brief Computes a maximal step size that is collision free.
///
/// All vertices in V0 and V1 will be considered for collisions, so V0 and
/// V1 should be only the surface vertices. The edges and face should be only
/// for the surface elements.
///
/// @note Assumes the trajectory is linear.
///
/// @param[in] mesh The surface of the contact mesh.
/// @param[in] V0
///     Vertex positions at start as rows of a matrix. Assumes V0 is
///     intersection free.
/// @param[in] V1 Surface vertex positions at end as rows of a matrix.
/// @returns A step-size \f$\in [0, 1]\f$ that is collision free. A value of 1.0
/// if a full step and 0.0 is no step.
double compute_collision_free_stepsize(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const BroadPhaseMethod& method = BroadPhaseMethod::HASH_GRID,
    const double tolerance = 1e-6,
    const long max_iterations = 1e7);

/// @brief Computes a maximal step size that is collision free using a set of
/// collision candidates.
///
/// @note Assumes the trajectory is linear.
///
/// @param[in] candidates Set of candidates to check for collisions.
/// @param[in] mesh The surface of the contact mesh.
/// @param[in] V0
///     Vertex positions at start as rows of a matrix. Assumes V0 is
///     intersection free.
/// @param[in] V1 Surface vertex positions at end as rows of a matrix.
////
/// @returns A step-size \f$\in [0, 1]\f$ that is collision free. A value of 1.0
/// if a full step and 0.0 is no step.
double compute_collision_free_stepsize(
    const Candidates& candidates,
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const double tolerance = 1e-6,
    const long max_iterations = 1e7);

///////////////////////////////////////////////////////////////////////////////
// Utilities

/// @brief Computes the minimum distance between any non-adjacent elements.
///
/// @param[in] mesh The surface of the contact mesh.
/// @param[in] V Surface vertex positions as rows of a matrix.
/// @returns The minimum distance between any non-adjacent elements.
double compute_minimum_distance(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const Constraints& constraint_set);

/// @brief Determine if the mesh has self intersections.
///
/// @param[in] mesh The surface of the contact mesh.
/// @param[in] V Surface vertex positions as rows of a matrix.
bool has_intersections(const CollisionMesh& mesh, const Eigen::MatrixXd& V);

} // namespace ipc
