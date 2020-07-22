#pragma once

#include <Eigen/Core>

#include <spatial_hash/collision_candidate.hpp>

namespace ipc {

void compute_constraint_set(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    double dhat_squared,
    ccd::Candidates& constraint_set);

double compute_barrier_potential(
    const Eigen::MatrixXd& V_rest,
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const ccd::Candidates& constraint_set,
    double dhat_squared,
    double barrier_stiffness);

bool is_step_collision_free(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F);

double compute_intersection_free_stepsize(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F);

double compute_minimum_distance(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F);

// double compute_friction_potential(
//     const Eigen::MatrixXd& V_rest,
//     const Eigen::MatrixXd& V,
//     const Eigen::MatrixXi& E,
//     const Eigen::MatrixXi& F,
//     const Candidates& constraint_set,
//     double dhat_squared,
//     double barrier_stiffness);

} // namespace ipc
