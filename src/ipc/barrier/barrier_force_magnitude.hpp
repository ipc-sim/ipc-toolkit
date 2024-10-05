#pragma once

#include <ipc/barrier/barrier.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

/// @brief Compute the magnitude of the force due to a barrier.
/// @param distance_squared The squared distance between elements.
/// @param barrier The barrier function.
/// @param dhat The activation distance of the barrier.
/// @param barrier_stiffness The stiffness of the barrier.
/// @param dmin The minimum distance offset to the barrier.
/// @return The magnitude of the force.
double barrier_force_magnitude(
    const double distance_squared,
    const Barrier& barrier,
    const double dhat,
    const double barrier_stiffness,
    const double dmin = 0);

/// @brief Compute the gradient of the magnitude of the force due to a barrier.
/// @param distance_squared The squared distance between elements.
/// @param distance_squared_gradient The gradient of the squared distance.
/// @param barrier The barrier function.
/// @param dhat The activation distance of the barrier.
/// @param barrier_stiffness The stiffness of the barrier.
/// @param dmin The minimum distance offset to the barrier.
/// @return The gradient of the force.
VectorMax12d barrier_force_magnitude_gradient(
    const double distance_squared,
    const VectorMax12d& distance_squared_gradient,
    const Barrier& barrier,
    const double dhat,
    const double barrier_stiffness,
    const double dmin = 0);

} // namespace ipc
