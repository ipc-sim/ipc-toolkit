#pragma once

#include <Eigen/Core>

#include <ipc/broad_phase/collision_candidate.hpp>

void brute_force_comparison(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const Eigen::VectorXi& group_ids,
    ipc::Candidates& candidates,
    const double inflation_radius);
