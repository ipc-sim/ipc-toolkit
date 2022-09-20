#include <catch2/catch.hpp>

#include <broad_phase/brute_force_comparison.hpp>

#include <tbb/parallel_sort.h>

#include <ipc/broad_phase/broad_phase.hpp>

#include <nlohmann/json.hpp>

#include <fstream>

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
        construct_collision_candidates(
            mesh, V0, V1, bf_candidates, inflation_radius,
            BroadPhaseMethod::BRUTE_FORCE);
        // save_candidates("bf_candidates.json", bf_candidates);
    }

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

void save_candidates(
    const std::string& filename, const ipc::Candidates& candidates)
{
    std::vector<std::array<long, 2>> ev_candidates;
    ev_candidates.reserve(candidates.ev_candidates.size());
    for (const auto& ev : candidates.ev_candidates) {
        ev_candidates.push_back({ { ev.edge_index, ev.vertex_index } });
    }

    std::vector<std::array<long, 2>> ee_candidates;
    ee_candidates.reserve(candidates.ee_candidates.size());
    for (const auto& ee : candidates.ee_candidates) {
        ee_candidates.push_back({ { ee.edge0_index, ee.edge1_index } });
    }

    std::vector<std::array<long, 2>> fv_candidates;
    fv_candidates.reserve(candidates.fv_candidates.size());
    for (const auto& fv : candidates.fv_candidates) {
        fv_candidates.push_back({ { fv.face_index, fv.vertex_index } });
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
