#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <ipc/collision_constraint.hpp>
// NOTE: Include this so the user can just include ipc.hpp
#include <ipc/friction/friction.hpp>

/// Incremental Potential Contact functions
namespace ipc {

/// @brief Construct a set of constraints used to compute the barrier potential.
///
/// @note V can either be the surface vertices or the entire mesh vertices.
/// The edges and face should be only for the surface elements.
///
/// @param[in] V Vertex positions as rows of a matrix.
/// @param[in] E Edges as rows of indicies into V.
/// @param[in] F Triangular faces as rows of indicies into V.
/// @param[in] dhat The activation distance of the barrier.
/// @param[out] constraint_set The constructed set of constraints.
/// @param[in] ignore_codimensional_vertices Ignores vertices not connected to
///                                          edges.
/// @param[in] vertex_group_ids A group ID per vertex such that vertices with
///                             the same group id do not collide. An empty
///                             vector implies all vertices can collide with all
///                             other vertices.
void construct_constraint_set(
    const Eigen::MatrixXd& V_rest,
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    double dhat,
    Constraints& constraint_set,
    bool ignore_codimensional_vertices = true,
    const Eigen::VectorXi& vertex_group_ids = Eigen::VectorXi(),
    const Eigen::MatrixXi& F2E = Eigen::MatrixXi(),
    double dmin = 0);

/// @brief Construct a set of constraints used to compute the barrier potential.
///
/// @note V can either be the surface vertices or the entire mesh vertices.
/// The edges and face should be only for the surface elements.
///
/// @param[in] candidates Distance candidates from which the constraint set is
///                       built.
/// @param[in] V Vertex positions as rows of a matrix.
/// @param[in] E Edges as rows of indicies into V.
/// @param[in] F Triangular faces as rows of indicies into V.
/// @param[in] dhat The activation distance of the barrier.
/// @param[out] constraint_set The constructed set of constraints.
void construct_constraint_set(
    const Candidates& candidates,
    const Eigen::MatrixXd& V_rest,
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    double dhat,
    Constraints& constraint_set,
    const Eigen::MatrixXi& F2E = Eigen::MatrixXi(),
    double dmin = 0);

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
    double dhat,
    double dmin = 0);

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
    double dhat,
    double dmin = 0);

/// @brief Compute the hessian of the barrier potential.
///
/// @param[in] V Vertex positions as rows of a matrix.
/// @param[in] E Edges as rows of indicies into V.
/// @param[in] F Triangular faces as rows of indicies into V.
/// @param[in] constraint_set The set of constraints.
/// @param[in] dhat The activation distance of the barrier.
/// @returns The hessian of all barrier potentials (not scaled by the barrier
/// stiffness).
Eigen::SparseMatrix<double> compute_barrier_potential_hessian(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const Constraints& constraint_set,
    double dhat,
    double dmin = 0,
    bool project_to_psd = true);

/// @brief Determine if the step is collision free.
///
/// @note Assumes the trajectory is linear.
/// @note V can either be the surface vertices or the entire mesh vertices.
/// The edges and face should be only for the surface elements.
///
/// @param[in] V0 Vertex positions at start as rows of a matrix.
/// @param[in] V1 Vertex positions at end as rows of a matrix.
/// @param[in] E Edges as rows of indicies into V.
/// @param[in] F Triangular faces as rows of indicies into V.
/// @param[in] ignore_codimensional_vertices Ignores vertices not connected to
///                                          edges.
/// @param[in] vertex_group_ids A group ID per vertex such that vertices with
///                             the same group id do not collide. An empty
///                             vector implies all vertices can collide with all
///                             other vertices.
/// @returns True if <b>any</b> collisions occur.
bool is_step_collision_free(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    bool ignore_codimensional_vertices = true,
    const Eigen::VectorXi& vertex_group_ids = Eigen::VectorXi());

/// @brief Computes a maximal step size that is collision free.
///
/// @note Assumes V0 is intersection free.
/// @note Assumes the trajectory is linear.
/// @note A value of 1.0 if a full step and 0.0 is no step.
///
/// @param[in] V0 Vertex positions at start as rows of a matrix.
/// @param[in] V1 Vertex positions at end as rows of a matrix.
/// @param[in] E Edges as rows of indicies into V.
/// @param[in] F Triangular faces as rows of indicies into V.
/// @param[in] ignore_codimensional_vertices Ignores vertices not connected to
///                                          edges.
/// @param[in] vertex_group_ids A group ID per vertex such that vertices with
///                             the same group id do not collide. An empty
///                             vector implies all vertices can collide with all
///                             other vertices.
/// @returns A step-size \f$\in [0, 1]\f$ that is collision free.
double compute_collision_free_stepsize(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    bool ignore_codimensional_vertices = true,
    const Eigen::VectorXi& vertex_group_ids = Eigen::VectorXi());

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

} // namespace ipc
