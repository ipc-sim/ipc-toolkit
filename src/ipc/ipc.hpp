#pragma once

// NOTE: Include this so the user can just include ipc.hpp
#include <ipc/collisions/constraints.hpp>
#include <ipc/friction/friction.hpp>
#include <ipc/broad_phase/broad_phase.hpp>
#include <ipc/collision_mesh.hpp>

#include <Eigen/Core>
#include <Eigen/Sparse>

/// Incremental Potential Contact functions
namespace ipc {

/// @brief Compute the barrier potential for a given constraint set.
/// @param[in] mesh The collision mesh.
/// @param[in] V Vertices of the collision mesh.
/// @param[in] constraint_set The set of constraints.
/// @param[in] dhat The activation distance of the barrier.
/// @returns The sum of all barrier potentials (not scaled by the barrier stiffness).
double compute_barrier_potential(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const Constraints& constraint_set,
    const double dhat);

/// @brief Compute the gradient of the barrier potential.
/// @param[in] mesh The collision mesh.
/// @param[in] V Vertices of the collision mesh.
/// @param[in] constraint_set The set of constraints.
/// @param[in] dhat The activation distance of the barrier.
/// @returns The gradient of all barrier potentials (not scaled by the barrier stiffness). This will have a size of |V|.
Eigen::VectorXd compute_barrier_potential_gradient(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const Constraints& constraint_set,
    const double dhat);

/// @brief Compute the hessian of the barrier potential.
/// @param[in] mesh The collision mesh.
/// @param[in] V Vertices of the collision mesh.
/// @param[in] constraint_set The set of constraints.
/// @param[in] dhat The activation distance of the barrier.
/// @param[in] project_hessian_to_psd Make sure the hessian is positive semi-definite.
/// @returns The hessian of all barrier potentials (not scaled by the barrier stiffness). This will have a size of |V|x|V|.
Eigen::SparseMatrix<double> compute_barrier_potential_hessian(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const Constraints& constraint_set,
    const double dhat,
    const bool project_hessian_to_psd = true);

/// @brief Compute the barrier shape derivative.
/// @param[in] mesh The collision mesh.
/// @param[in] V Vertices of the collision mesh.
/// @param[in] constraint_set The set of constraints.
/// @param[in] dhat The activation distance of the barrier.
/// @returns The derivative of the force with respect to X, the rest positions.
Eigen::SparseMatrix<double> compute_barrier_shape_derivative(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const Constraints& constraint_set,
    const double dhat);

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

/// @brief Computes the minimum distance between any non-adjacent elements.
/// @param[in] mesh The collision mesh.
/// @param[in] V Vertices of the collision mesh.
/// @returns The minimum distance between any non-adjacent elements.
double compute_minimum_distance(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const Constraints& constraint_set);

/// @brief Determine if the mesh has self intersections.
/// @param[in] mesh The collision mesh.
/// @param[in] V Vertices of the collision mesh.
/// @return A boolean for if the mesh has intersections.
bool has_intersections(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const BroadPhaseMethod method = BroadPhaseMethod::HASH_GRID);

} // namespace ipc
