// Fucntions for computing the initial and updated barrier stiffnesses.

#pragma once

#include <Eigen/Core>

#include <ipc/collision_constraint.hpp>
#include <ipc/collision_mesh.hpp>

namespace ipc {

/// Compute an inital barrier stiffness using the barrier potential gradient.
double initial_barrier_stiffness(
    double bbox_diagonal,
    double dhat,
    double average_mass,
    const Eigen::VectorXd& grad_energy,
    const Eigen::VectorXd& grad_barrier,
    double& max_barrier_stiffness,
    double min_barrier_stiffness_scale = 1e11,
    double dmin = 0);

/// Update the barrier stiffness if the distance is decreasing and less than
/// dhat_epsilon_scale * diag.
void update_barrier_stiffness(
    double prev_min_distance,
    double min_distance,
    double max_barrier_stiffness,
    double& barrier_stiffness,
    double bbox_diagonal,
    double dhat_epsilon_scale = 1e-9,
    double dmin = 0);

/// Update the barrier stiffness if the distance is decreasing and less than
/// dhat_epsilon_scale * diag.
void update_barrier_stiffness(
    const Constraints& constraint_set,
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    double prev_min_distance,
    double& min_distance,
    double max_barrier_stiffness,
    double& barrier_stiffness,
    double dhat_epsilon_scale = 1e-9,
    double dmin = 0);

} // namespace ipc
