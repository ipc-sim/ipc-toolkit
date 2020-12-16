#pragma once
#include <ipc/friction/friction.hpp>

namespace ipc {

template <typename T>
T compute_friction_potential(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixX<T>& V1,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const FrictionConstraints& friction_constraint_set,
    double epsv_times_h)
{
    // absolute linear dislacement of each point
    Eigen::MatrixX<T> U = V1 - V0.cast<T>();

    T potential(0);

    for (const auto& vv_constraint : friction_constraint_set.vv_constraints) {
        potential += vv_constraint.compute_potential(U, E, F, epsv_times_h);
    }

    for (const auto& ev_constraint : friction_constraint_set.ev_constraints) {
        potential += ev_constraint.compute_potential(U, E, F, epsv_times_h);
    }

    for (const auto& ee_constraint : friction_constraint_set.ee_constraints) {
        potential += ee_constraint.compute_potential(U, E, F, epsv_times_h);
    }

    for (const auto& fv_constraint : friction_constraint_set.fv_constraints) {
        potential += fv_constraint.compute_potential(U, E, F, epsv_times_h);
    }

    return potential;
}

} // namespace ipc
