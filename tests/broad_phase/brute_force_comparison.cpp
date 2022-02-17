#include <catch2/catch.hpp>

#include <broad_phase/brute_force_comparison.hpp>

#include <tbb/parallel_sort.h>

#include <ipc/broad_phase/broad_phase.hpp>

void brute_force_comparison(
    const ipc::CollisionMesh& mesh,
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    ipc::Candidates& candidates,
    const double inflation_radius)
{
    using namespace ipc;

    Candidates bf_candidates;
    construct_collision_candidates(
        mesh, V0, V1, bf_candidates, inflation_radius,
        BroadPhaseMethod::BRUTE_FORCE);

    CAPTURE(candidates.size(), bf_candidates.size());

    CAPTURE("EV");
    brute_force_comparison(
        mesh, V0, V1, candidates.ev_candidates, bf_candidates.ev_candidates);
    CAPTURE("EE");
    brute_force_comparison(
        mesh, V0, V1, candidates.ee_candidates, bf_candidates.ee_candidates);
    CAPTURE("FV");
    brute_force_comparison(
        mesh, V0, V1, candidates.fv_candidates, bf_candidates.fv_candidates);
}

template <typename Candidate>
void brute_force_comparison(
    const ipc::CollisionMesh& mesh,
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    std::vector<Candidate>& candidates,
    std::vector<Candidate>& bf_candidates)
{
    CHECK(candidates.size() <= bf_candidates.size());

    tbb::parallel_sort(candidates.begin(), candidates.end());
    tbb::parallel_sort(bf_candidates.begin(), bf_candidates.end());

    int ci = 0;
    for (int bf_ci = 0; bf_ci < bf_candidates.size(); bf_ci++) {
        if (candidates.size() <= ci || bf_candidates[bf_ci] != candidates[ci]) {
            // Perform CCD to make sure the candidate is not a collision
            double toi;
            bool hit = bf_candidates[bf_ci].ccd(
                V0, V1, mesh.edges(), mesh.faces(), toi, /*tmax=*/1.0,
                /*tolerance=*/1e-6, /*max_iterations=*/1e7,
                /*conservative_rescaling=*/1.0);
            CHECK(!hit); // Check for FN
        } else {
            ci++;
        }
    }

    CHECK(ci >= candidates.size());
}

template void brute_force_comparison<ipc::EdgeVertexCandidate>(
    const ipc::CollisionMesh& mesh,
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    std::vector<ipc::EdgeVertexCandidate>& candidates,
    std::vector<ipc::EdgeVertexCandidate>& bf_candidates);
template void brute_force_comparison<ipc::EdgeEdgeCandidate>(
    const ipc::CollisionMesh& mesh,
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    std::vector<ipc::EdgeEdgeCandidate>& candidates,
    std::vector<ipc::EdgeEdgeCandidate>& bf_candidates);
template void brute_force_comparison<ipc::FaceVertexCandidate>(
    const ipc::CollisionMesh& mesh,
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    std::vector<ipc::FaceVertexCandidate>& candidates,
    std::vector<ipc::FaceVertexCandidate>& bf_candidates);