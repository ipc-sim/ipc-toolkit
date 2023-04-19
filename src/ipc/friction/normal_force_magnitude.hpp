#pragma once

#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

double compute_normal_force_magnitude(
    double distance_squared,
    double dhat,
    double barrier_stiffness,
    double dmin = 0);

VectorMax12d compute_normal_force_magnitude_gradient(
    double distance_squared,
    const Eigen::VectorXd& distance_squared_gradient,
    double dhat,
    double barrier_stiffness,
    double dmin = 0);

} // namespace ipc
