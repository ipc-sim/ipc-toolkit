#include "brute_force_comparison.hpp"

#include <catch2/catch_test_macros.hpp>

#include <tbb/parallel_sort.h>
#include <nlohmann/json.hpp>

void save_candidates(
    const std::string& filename, const ipc::Candidates& candidates);
bool load_candidates(const std::string& filename, ipc::Candidates& candidates);

void brute_force_comparison(
    const ipc::CollisionMesh& mesh,
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    ipc::Candidates& candidates,
    const double inflation_radius,
    const std::string& cached_bf_candidates)
{
    using namespace ipc;

    Candidates bf_candidates;
    if (cached_bf_candidates.empty()
        || !load_candidates(cached_bf_candidates, bf_candidates)) {
        bf_candidates.build(
            mesh, V0, V1, inflation_radius, BroadPhaseMethod::BRUTE_FORCE);
        if (!cached_bf_candidates.empty()) {
            save_candidates(cached_bf_candidates, bf_candidates);
        }
    }

    CAPTURE(candidates.size(), bf_candidates.size());

    {
        INFO("EV");
        brute_force_comparison(
            mesh, V0, V1, candidates.ev_candidates,
            bf_candidates.ev_candidates);
    }

    {
        INFO("EE");
        brute_force_comparison(
            mesh, V0, V1, candidates.ee_candidates,
            bf_candidates.ee_candidates);
    }

    {
        INFO("FV");
        brute_force_comparison(
            mesh, V0, V1, candidates.fv_candidates,
            bf_candidates.fv_candidates);
    }
}

void save_candidates(
    const std::string& filename, const ipc::Candidates& candidates)
{
    std::vector<std::array<long, 2>> ev_candidates;
    ev_candidates.reserve(candidates.ev_candidates.size());
    for (const auto& ev : candidates.ev_candidates) {
        ev_candidates.push_back({ { ev.edge_id, ev.vertex_id } });
    }

    std::vector<std::array<long, 2>> ee_candidates;
    ee_candidates.reserve(candidates.ee_candidates.size());
    for (const auto& ee : candidates.ee_candidates) {
        ee_candidates.push_back({ { ee.edge0_id, ee.edge1_id } });
    }

    std::vector<std::array<long, 2>> fv_candidates;
    fv_candidates.reserve(candidates.fv_candidates.size());
    for (const auto& fv : candidates.fv_candidates) {
        fv_candidates.push_back({ { fv.face_id, fv.vertex_id } });
    }

    nlohmann::json out;
    out["ev_candidates"] = ev_candidates;
    out["ee_candidates"] = ee_candidates;
    out["fv_candidates"] = fv_candidates;

    std::ofstream f(filename);
    f << out;
}

bool load_candidates(const std::string& filename, ipc::Candidates& candidates)
{
    std::ifstream f(filename);
    if (!f)
        return false;

    nlohmann::json in;
    f >> in;

    std::vector<std::array<long, 2>> ev_candidates = in["ev_candidates"];
    candidates.ev_candidates.reserve(ev_candidates.size());
    for (const auto& [ei, vi] : ev_candidates) {
        candidates.ev_candidates.emplace_back(ei, vi);
    }

    std::vector<std::array<long, 2>> ee_candidates = in["ee_candidates"];
    candidates.ee_candidates.reserve(ee_candidates.size());
    for (const auto& [e0i, e1i] : ee_candidates) {
        candidates.ee_candidates.emplace_back(e0i, e1i);
    }

    std::vector<std::array<long, 2>> fv_candidates = in["fv_candidates"];
    fv_candidates.reserve(candidates.fv_candidates.size());
    for (const auto& [fi, vi] : fv_candidates) {
        candidates.fv_candidates.emplace_back(fi, vi);
    }

    return true;
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

    for (const Candidate& bf_candidate : bf_candidates) {
        CAPTURE(bf_candidate);
        const bool found = std::binary_search(
            candidates.begin(), candidates.end(), bf_candidate);
        if (!found) {
            // Perform CCD to make sure the candidate is not a collision
            double toi;
            const bool hit = bf_candidate.ccd(
                bf_candidate.dof(V0, mesh.edges(), mesh.faces()),
                bf_candidate.dof(V1, mesh.edges(), mesh.faces()), //
                toi, /*min_distance=*/0, /*tmax=*/1.0,
                ipc::DEFAULT_CCD_TOLERANCE, ipc::DEFAULT_CCD_MAX_ITERATIONS,
                /*conservative_rescaling=*/1.0);
            CHECK(!hit); // Check for FN
        }
    }
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
