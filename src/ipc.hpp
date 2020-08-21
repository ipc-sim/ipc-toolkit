#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <spatial_hash/collision_candidate.hpp>

namespace ipc {

/// @brief Construct a set of constraints used to compute the barrier potential.
///
/// @note V can either be the surface vertices or the entire mesh vertices.
/// The edges and face should be only for the surface elements.
///
/// @param[in] V Vertex positions as rows of a matrix.
/// @param[in] E Edges as rows of indicies into V.
/// @param[in] F Triangular faces as rows of indicies into V.
/// @param[in] dhat_squared The activation distance squared of the barrier.
/// @param[out] constraint_set The constructed set of constraints.
/// @param[in] ignore_internal_vertices Ignores vertices not connected to edges E
void construct_constraint_set(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    double dhat_squared,
    ccd::Candidates& constraint_set,
    bool ignore_internal_vertices = true);

/// @brief Construct a set of constraints used to compute the barrier potential.
///
/// @note V can either be the surface vertices or the entire mesh vertices.
/// The edges and face should be only for the surface elements.
///
/// @param[in] V_rest Vertex positions at rest as rows of a matrix.
/// @param[in] V Vertex positions as rows of a matrix.
/// @param[in] E Edges as rows of indicies into V.
/// @param[in] F Triangular faces as rows of indicies into V.
/// @param[in] constraint_set The set of constraints.
/// @param[in] dhat_squared The activation distance squared of the barrier.
/// @returns The sum of all barrier potentials.
double compute_barrier_potential(
    const Eigen::MatrixXd& V_rest,
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const ccd::Candidates& constraint_set,
    double dhat_squared);

Eigen::VectorXd compute_barrier_potential_gradient(
    const Eigen::MatrixXd& V_rest,
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const ccd::Candidates& constraint_set,
    double dhat_squared);

Eigen::SparseMatrix<double> compute_barrier_potential_hessian(
    const Eigen::MatrixXd& V_rest,
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const ccd::Candidates& constraint_set,
    double dhat_squared);

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
/// @returns True if <b>any</b> collisions occur.
bool is_step_collision_free(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F);

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
/// @returns A step-size \f$\in [0, 1]\f$ that is collision free.
double compute_collision_free_stepsize(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F);

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
    const ccd::Candidates& constraint_set);

// double compute_friction_potential(
//     const Eigen::MatrixXd& V_prev,
//     const Eigen::MatrixXd& V,
//     const Eigen::MatrixXi& E,
//     const Eigen::MatrixXi& F,
//     const Candidates& constraint_set,
//     double epsilon_v);

} // namespace ipc
