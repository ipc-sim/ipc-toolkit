#pragma once

#include <ipc/candidates/candidate_vector.hpp>
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
    ipc::CandidateVector<Candidate>& candidates,
    ipc::CandidateVector<Candidate>& bf_candidates);
