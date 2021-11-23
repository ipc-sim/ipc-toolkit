#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <ipc/collision_constraint.hpp>
// NOTE: Include this so the user can just include ipc.hpp
#include <ipc/friction/friction.hpp>
#include <ipc/broad_phase/broad_phase.hpp>

/// Incremental Potential Contact functions
namespace ipc {

/// @brief Construct a set of constraints used to compute the barrier potential.
///
/// All vertices in V will be considered for collisions, so V should be only the
/// surface vertices. The edges and face should be only for the surface
/// elements.
///
/// @param[in] V Vertex positions as rows of a matrix.
/// @param[in] E Edges as rows of indicies into V.
/// @param[in] F Triangular faces as rows of indicies into V.
/// @param[in] dhat The activation distance of the barrier.
/// @param[out] constraint_set
///     The constructed set of constraints (any existing constraints will be
///     cleared).
/// @param[in] F2E Map from F edges to rows of E.
/// @param[in] dmin Minimum distance.
/// @param[in] method Broad-phase method to use.
/// @param[in] can_collide
///     A function that takes two vertex IDs (row numbers in V) and returns true
///     if the vertices (and faces or edges containing the vertices) can
///     collide. By default all primitives can collide with all other
///     primitives.
void construct_constraint_set(
    const Eigen::MatrixXd& V_rest,
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const double dhat,
    Constraints& constraint_set,
    const Eigen::MatrixXi& F2E = Eigen::MatrixXi(),
    const double dmin = 0,
    const BroadPhaseMethod& method = BroadPhaseMethod::HASH_GRID,
    const std::function<bool(size_t, size_t)>& can_collide =
        [](size_t, size_t) { return true; });

/// @brief Construct a set of constraints used to compute the barrier potential.
///
/// V can either be the surface vertices or the entire mesh vertices. The edges
/// and face should be only for the surface elements.
///
/// @param[in] V Vertex positions as rows of a matrix.
/// @param[in] codim_V
///     Vertex indices of codimensional points. Pass `Eigen::VectorXi()` to
///     ignore all codimensional vertices (vertices not connected to an edge).
/// @param[in] E Edges as rows of indicies into V.
/// @param[in] F Triangular faces as rows of indicies into V.
/// @param[in] dhat The activation distance of the barrier.
/// @param[out] constraint_set
///     The constructed set of constraints (any existing constraints will be
///     cleared).
/// @param[in] F2E Map from F edges to rows of E.
/// @param[in] dmin Minimum distance.
/// @param[in] method Broad-phase method to use.
/// @param[in] can_collide
///     A function that takes two vertex IDs (row numbers in V) and returns true
///     if the vertices (and faces or edges containing the vertices) can
///     collide. By default all primitives can collide with all other
///     primitives.
void construct_constraint_set(
    const Eigen::MatrixXd& V_rest,
    const Eigen::MatrixXd& V,
    const Eigen::VectorXi& codim_V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const double dhat,
    Constraints& constraint_set,
    const Eigen::MatrixXi& F2E = Eigen::MatrixXi(),
    const double dmin = 0,
    const BroadPhaseMethod& method = BroadPhaseMethod::HASH_GRID,
    const std::function<bool(size_t, size_t)>& can_collide =
        [](size_t, size_t) { return true; });

/// @brief Construct a set of constraints used to compute the barrier potential.
///
/// V can either be the surface vertices or the entire mesh vertices. The edges
/// and face should be only for the surface elements.
///
/// @param[in]  candidates  Distance candidates from which the constraint set is
///                         built.
/// @param[in]  V  Vertex positions as rows of a matrix.
/// @param[in]  E  Edges as rows of indicies into V.
/// @param[in]  F  Triangular faces as rows of indicies into V.
/// @param[in]  dhat  The activation distance of the barrier.
/// @param[out] constraint_set
///     The constructed set of constraints (any existing constraints will be
///     cleared).
/// @param[in]  F2E  Map from F edges to rows of E.
/// @param[in]  dmin  Minimum distance.
void construct_constraint_set(
    const Candidates& candidates,
    const Eigen::MatrixXd& V_rest,
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const double dhat,
    Constraints& constraint_set,
    const Eigen::MatrixXi& F2E = Eigen::MatrixXi(),
    const double dmin = 0);

/// @brief Compute the barrier potential for a given constraint set.
///
/// @param[in] V Vertex positions as rows of a matrix.
/// @param[in] E Edges as rows of indicies into V.
/// @param[in] F Triangular faces as rows of indicies into V.
/// @param[in] constraint_set The set of constraints.
/// @param[in] dhat The activation distance of the barrier.
/// @returns The sum of all barrier potentials (not scaled by the barrier
/// stiffness).
double compute_barrier_potential(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const Constraints& constraint_set,
    const double dhat);

/// @brief Compute the gradient of the barrier potential.
///
/// @param[in] V Vertex positions as rows of a matrix.
/// @param[in] E Edges as rows of indicies into V.
/// @param[in] F Triangular faces as rows of indicies into V.
/// @param[in] constraint_set The set of constraints.
/// @param[in] dhat The activation distance of the barrier.
/// @returns The gradient of all barrier potentials (not scaled by the barrier
/// stiffness).
Eigen::VectorXd compute_barrier_potential_gradient(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const Constraints& constraint_set,
    const double dhat);

/// @brief Compute the hessian of the barrier potential.
///
/// @param[in] V Vertex positions as rows of a matrix.
/// @param[in] E Edges as rows of indicies into V.
/// @param[in] F Triangular faces as rows of indicies into V.
/// @param[in] constraint_set The set of constraints.
/// @param[in] dhat The activation distance of the barrier.
/// @param[in] project_hessian_to_psd Make sure the hessian is positive
///            semi-definite.
/// @returns The hessian of all barrier potentials (not scaled by the barrier
/// stiffness).
Eigen::SparseMatrix<double> compute_barrier_potential_hessian(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const Constraints& constraint_set,
    const double dhat,
    const bool project_hessian_to_psd = true);

///////////////////////////////////////////////////////////////////////////////
// Collision detection

/// @brief Determine if the step is collision free.
///
/// All vertices in V0 and V1 will be considered for collisions, so V0 and
/// V1 should be only the surface vertices. The edges and face should be only
/// for the surface elements.
///
/// @note Assumes the trajectory is linear.
///
/// @param[in] V0 Vertex positions at start as rows of a matrix.
/// @param[in] V1 Vertex positions at end as rows of a matrix.
/// @param[in] E Edges as rows of indicies into V.
/// @param[in] F Triangular faces as rows of indicies into V.
/// @param[in] can_collide
///     A function that takes two vertex IDs (row numbers in F) and returns true
///     if the vertices (and faces or edges containing the vertices) can
///     collide. By default all primitives can collide with all other primiti
/// @returns True if <b>any</b> collisions occur.
bool is_step_collision_free(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const BroadPhaseMethod& method = BroadPhaseMethod::HASH_GRID,
    const double tolerance = 1e-6,
    const long max_iterations = 1e7,
    const std::function<bool(size_t, size_t)>& can_collide =
        [](size_t, size_t) { return true; });

/// @brief Determine if the step is collision free.
///
/// V can either be the surface vertices or the entire mesh vertices. The edges
/// and face should be only for the surface elements.
///
/// @note Assumes the trajectory is linear.
///
/// @param[in] V0 Vertex positions at start as rows of a matrix.
/// @param[in] V1 Vertex positions at end as rows of a matrix.
/// @param[in] codim_V Vertex indices of codimensional points. Pass
///     Eigen::VectorXi() to ignore all codimensional vertices (vertices not
///     connected to an edge).
/// @param[in] E Edges as rows of indicies into V.
/// @param[in] F Triangular faces as rows of indicies into V.
/// @param[in] can_collide
///     A function that takes two vertex IDs (row numbers in F) and returns true
///     if the vertices (and faces or edges containing the vertices) can
///     collide. By default all primitives can collide with all other primiti
/// @returns True if <b>any</b> collisions occur.
bool is_step_collision_free(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::VectorXi& codim_V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const BroadPhaseMethod& method = BroadPhaseMethod::HASH_GRID,
    const double tolerance = 1e-6,
    const long max_iterations = 1e7,
    const std::function<bool(size_t, size_t)>& can_collide =
        [](size_t, size_t) { return true; });

/// @brief Determine if the step is collision free from a set of candidates.
///
/// V can either be the surface vertices or the entire mesh vertices. The edges
/// and face should be only for the surface elements.
///
/// @note Assumes the trajectory is linear.
///
/// @param[in] candidates Set of candidates to check for collisions.
/// @param[in] V0 Vertex positions at start as rows of a matrix.
/// @param[in] V1 Vertex positions at end as rows of a matrix.
/// @param[in] E Edges as rows of indicies into V.
/// @param[in] F Triangular faces as rows of indicies into V.
///
/// @returns True if <b>any</b> collisions occur.
bool is_step_collision_free(
    const Candidates& candidates,
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
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
/// @param[in] V0
///     Vertex positions at start as rows of a matrix. Assumes V0 is
///     intersection free.
/// @param[in] V1 Vertex positions at end as rows of a matrix.
/// @param[in] E Edges as rows of indicies into V.
/// @param[in] F Triangular faces as rows of indicies into V.
/// @param[in] can_collide
///     A function that takes two vertex IDs (row numbers in F) and returns true
///     if the vertices (and faces or edges containing the vertices) can
///     collide. By default all primitives can collide with all other
///     primitives.
/// @returns A step-size \f$\in [0, 1]\f$ that is collision free. A value of 1.0
/// if a full step and 0.0 is no step.
double compute_collision_free_stepsize(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const BroadPhaseMethod& method = BroadPhaseMethod::HASH_GRID,
    const double tolerance = 1e-6,
    const long max_iterations = 1e7,
    const std::function<bool(size_t, size_t)>& can_collide =
        [](size_t, size_t) { return true; });

/// @brief Computes a maximal step size that is collision free.
///
/// @note Assumes the trajectory is linear.
///
/// @param[in] V0
///     Vertex positions at start as rows of a matrix. Assumes V0 is
///     intersection free.
/// @param[in] V1 Vertex positions at end as rows of a matrix.
/// @param[in] codim_V Vertex indices of codimensional points. Pass
///     Eigen::VectorXi() to ignore all codimensional vertices (vertices not
///     connected to an edge).
/// @param[in] E Edges as rows of indicies into V.
/// @param[in] F Triangular faces as rows of indicies into V.
/// @param[in] can_collide
///     A function that takes two vertex IDs (row numbers in F) and returns true
///     if the vertices (and faces or edges containing the vertices) can
///     collide. By default all primitives can collide with all other
///     primitives.
/// @returns A step-size \f$\in [0, 1]\f$ that is collision free. A value of 1.0
/// if a full step and 0.0 is no step.
double compute_collision_free_stepsize(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::VectorXi& codim_V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const BroadPhaseMethod& method = BroadPhaseMethod::HASH_GRID,
    const double tolerance = 1e-6,
    const long max_iterations = 1e7,
    const std::function<bool(size_t, size_t)>& can_collide =
        [](size_t, size_t) { return true; });

/// @brief Computes a maximal step size that is collision free using a set of
/// collision candidates.
///
/// @note Assumes the trajectory is linear.
///
/// @param[in] candidates Set of candidates to check for collisions.
/// @param[in] V0
///     Vertex positions at start as rows of a matrix. Assumes V0 is
///     intersection free.
/// @param[in] V1 Vertex positions at end as rows of a matrix.
/// @param[in] E Edges as rows of indicies into V.
/// @param[in] F Triangular faces as rows of indicies into V.
////
/// @returns A step-size \f$\in [0, 1]\f$ that is collision free. A value of 1.0
/// if a full step and 0.0 is no step.
double compute_collision_free_stepsize(
    const Candidates& candidates,
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const double tolerance = 1e-6,
    const long max_iterations = 1e7);

///////////////////////////////////////////////////////////////////////////////
// Utilities

/// @brief Computes the minimum distance between any non-adjacent elements.
///
/// @param[in] V Vertex positions as rows of a matrix.
/// @param[in] E Edges as rows of indicies into V.
/// @param[in] F Triangular faces as rows of indicies into V.
/// @returns The minimum distance between any non-adjacent elements.
double compute_minimum_distance(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const Constraints& constraint_set);

/// @brief Determine if the mesh has self intersections.
///
/// @param[in] V Vertex positions as rows of a matrix.
/// @param[in] E Edges as rows of indicies into V.
/// @param[in] F Triangular faces as rows of indicies into V.
/// @param[in] can_collide
///     A function that takes two vertex IDs (row numbers in F) and returns true
///     if the vertices (and faces or edges containing the vertices) can
///     collide. By default all primitives can collide with all other
///     primitives.
bool has_intersections(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const std::function<bool(size_t, size_t)>& can_collide =
        [](size_t, size_t) { return true; });

} // namespace ipc
