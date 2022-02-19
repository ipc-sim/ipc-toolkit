#pragma once

#include <Eigen/Core>

#include <ipc/collision_mesh.hpp>
#include <ipc/broad_phase/collision_candidate.hpp>

void brute_force_comparison(
    const ipc::CollisionMesh& mesh,
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    ipc::Candidates& candidates,
    const double inflation_radius);

template <typename Candidate>
void brute_force_comparison(
    const ipc::CollisionMesh& mesh,
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    std::vector<Candidate>& candidates,
    std::vector<Candidate>& bf_candidates);