#pragma once

#include <ipc/barrier/barrier.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

double compute_normal_force_magnitude(
    const Barrier& barrier,
    const double distance_squared,
    const double dhat,
    const double barrier_stiffness,
    const double dmin = 0);

VectorMax12d compute_normal_force_magnitude_gradient(
    const Barrier& barrier,
    const double distance_squared,
    const Eigen::VectorXd& distance_squared_gradient,
    const double dhat,
    const double barrier_stiffness,
    const double dmin = 0);

} // namespace ipc
