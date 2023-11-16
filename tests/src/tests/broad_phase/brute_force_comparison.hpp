#pragma once

#include <ipc/collision_mesh.hpp>
#include <ipc/candidates/candidates.hpp>

void brute_force_comparison(
    const ipc::CollisionMesh& mesh,
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    ipc::Candidates& candidates,
    const double inflation_radius,
    const std::string& cached_bf_candidates = "");

template <typename Candidate>
void brute_force_comparison(
    const ipc::CollisionMesh& mesh,
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    std::vector<Candidate>& candidates,
    std::vector<Candidate>& bf_candidates);
