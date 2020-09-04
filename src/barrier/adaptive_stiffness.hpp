// Fucntions for computing the initial and updated barrier stiffnesses.

#pragma once

#include <Eigen/Core>

namespace ipc {

double intial_barrier_stiffness(
    const Eigen::MatrixXd& V_rest,
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    double dhat,
    double average_mass,
    const Eigen::VectorXd& grad_energy,
    double& max_barrier_stiffness,
    double min_barrier_stiffness_scale = 1e11);

void update_barrier_stiffness(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    double dhat,
    double prev_min_distance,
    double& min_distance,
    double max_barrier_stiffness,
    double& barrier_stiffness,
    double dhat_epsilon_scale);

} // namespace ipc
