// Fucntions for computing the initial and updated barrier stiffnesses.

#pragma once

#include <ipc/barrier/barrier.hpp>

#include <Eigen/Core>

namespace ipc {

/// @brief Compute an inital barrier stiffness using the barrier potential gradient.
/// @param[in] bbox_diagonal Length of the diagonal of the bounding box of the scene.
/// @param[in] barrier Barrier function.
/// @param[in] dhat Activation distance of the barrier.
/// @param[in] average_mass Average mass of all bodies.
/// @param[in] grad_energy Gradient of the elasticity energy function.
/// @param[in] grad_barrier Gradient of the barrier potential.
/// @param[out] max_barrier_stiffness Maximum stiffness of the barrier.
/// @param[in] min_barrier_stiffness_scale Scale used to premultiply the minimum barrier stiffness.
/// @param[in] dmin Minimum distance between elements.
/// @return The initial barrier stiffness.
double initial_barrier_stiffness(
    const double bbox_diagonal,
    const Barrier& barrier,
    const double dhat,
    const double average_mass,
    const Eigen::VectorXd& grad_energy,
    const Eigen::VectorXd& grad_barrier,
    double& max_barrier_stiffness,
    const double min_barrier_stiffness_scale = 1e11,
    const double dmin = 0);

/// @brief Update the barrier stiffness if the distance is decreasing and less than dhat_epsilon_scale * diag.
/// @param[in] prev_min_distance Previous minimum distance between elements.
/// @param[in] min_distance Current minimum distance between elements.
/// @param[in] max_barrier_stiffness Maximum stiffness of the barrier.
/// @param[in] barrier_stiffness Current barrier stiffness.
/// @param[in] bbox_diagonal Length of the diagonal of the bounding box of the scene.
/// @param[in] dhat_epsilon_scale Update if distance is less than this fraction of the diagonal.
/// @param[in] dmin Minimum distance between elements.
/// @return The updated barrier stiffness.
double update_barrier_stiffness(
    const double prev_min_distance,
    const double min_distance,
    const double max_barrier_stiffness,
    const double barrier_stiffness,
    const double bbox_diagonal,
    const double dhat_epsilon_scale = 1e-9,
    const double dmin = 0);

} // namespace ipc
