#include <tests/config.hpp>
#include <tests/utils.hpp>

#include <catch2/catch_test_macros.hpp>

#include <ipc/broad_phase/sweep_and_tiniest_queue.hpp>

using namespace ipc;

TEST_CASE("STQ Missing Features", "[broad_phase][stq]")
{
    Eigen::MatrixXd V0, V1;
    Eigen::MatrixXi E, F;

    REQUIRE(tests::load_mesh("cloth_ball92.ply", V0, E, F));
    REQUIRE(tests::load_mesh("cloth_ball93.ply", V1, E, F));

    double inflation_radius = 0;

#ifdef IPC_TOOLKIT_WITH_CUDA
    const BroadPhaseMethod method = GENERATE(
        BroadPhaseMethod::SWEEP_AND_TINIEST_QUEUE,
        BroadPhaseMethod::SWEEP_AND_TINIEST_QUEUE_GPU);
#else
    const BroadPhaseMethod method = BroadPhaseMethod::SWEEP_AND_TINIEST_QUEUE;
#endif

    std::shared_ptr<BroadPhase> stq = BroadPhase::make_broad_phase(method);
    stq->build(V0, V1, E, F, inflation_radius);

    Candidates candidates;
    stq->detect_collision_candidates(V0.cols(), candidates);

    CHECK(candidates.size() == 6'852'873);
    CHECK(candidates.vv_candidates.size() == 0);
    CHECK(candidates.ev_candidates.size() == 0);
    CHECK(candidates.ee_candidates.size() == 5'197'332);
    CHECK(candidates.fv_candidates.size() == 1'655'541);

    // Test missing features for code coverage

    try {
        std::vector<VertexVertexCandidate> vv_candidates;
        stq->detect_vertex_vertex_candidates(vv_candidates);
        FAIL("Should have thrown");
    } catch (const std::runtime_error& e) {
        SUCCEED(e.what());
    }

    try {
        std::vector<EdgeVertexCandidate> ev_candidates;
        stq->detect_edge_vertex_candidates(ev_candidates);
        FAIL("Should have thrown");
    } catch (const std::runtime_error& e) {
        SUCCEED(e.what());
    }

    try {
        std::vector<EdgeFaceCandidate> ef_candidates;
        stq->detect_edge_face_candidates(ef_candidates);
        FAIL("Should have thrown");
    } catch (const std::runtime_error& e) {
        SUCCEED(e.what());
    }

    try {
        std::vector<FaceFaceCandidate> ff_candidates;
        stq->detect_face_face_candidates(ff_candidates);
        FAIL("Should have thrown");
    } catch (const std::runtime_error& e) {
        SUCCEED(e.what());
    }

    stq->clear();
}
