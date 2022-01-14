#pragma once

#include <Eigen/Core>

namespace ipc {

double compute_normal_force_magnitude(
    double distance_squared,
    double dhat,
    double barrier_stiffness,
    double dmin = 0);


} // namespace ipc
