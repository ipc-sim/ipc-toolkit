#include <tests/config.hpp>
#include <tests/utils.hpp>

#include <catch2/catch_test_macros.hpp>

#include <ipc/broad_phase/sweep_and_prune.hpp>
#include <ipc/broad_phase/sweep_and_tiniest_queue.hpp>
#include <ipc/broad_phase/bvh.hpp>

using namespace ipc;

TEST_CASE("STQ All Cases", "[broad_phase][stq][cuda]")
{
#ifndef NDEBUG
    SKIP("'STQ All Cases' test is skipped in debug mode");
#endif

    Eigen::MatrixXd V0, V1;
    Eigen::MatrixXi E, F;

    REQUIRE(tests::load_mesh("cloth_ball92.ply", V0, E, F));
    REQUIRE(tests::load_mesh("cloth_ball93.ply", V1, E, F));

    double inflation_radius = 0;

#ifdef IPC_TOOLKIT_WITH_CUDA
    const std::shared_ptr<BroadPhase> broad_phase = GENERATE(
        std::static_pointer_cast<BroadPhase>(std::make_shared<SweepAndPrune>()),
        std::static_pointer_cast<BroadPhase>(
            std::make_shared<SweepAndTiniestQueue>()));
#else
    const std::shared_ptr<BroadPhase> broad_phase =
        std::make_shared<SweepAndPrune>();
#endif

    broad_phase->build(V0, V1, E, F, inflation_radius);

    Candidates candidates;
    broad_phase->detect_collision_candidates(V0.cols(), candidates);

    CHECK(candidates.size() == 6'852'873);
    CHECK(candidates.vv_candidates.size() == 0);
    CHECK(candidates.ev_candidates.size() == 0);
    CHECK(candidates.ee_candidates.size() == 5'197'332);
    CHECK(candidates.fv_candidates.size() == 1'655'541);

    std::vector<VertexVertexCandidate> vv_candidates;
    broad_phase->detect_vertex_vertex_candidates(vv_candidates);
    CHECK(vv_candidates.size() == 84'912);

    std::vector<EdgeVertexCandidate> ev_candidates;
    broad_phase->detect_edge_vertex_candidates(ev_candidates);
    CHECK(ev_candidates.size() == 1'666'926);

    std::vector<EdgeFaceCandidate> ef_candidates;
    broad_phase->detect_edge_face_candidates(ef_candidates);
    CHECK(ef_candidates.size() == 9'248'220);

    std::vector<FaceFaceCandidate> ff_candidates;
    broad_phase->detect_face_face_candidates(ff_candidates);
    CHECK(ff_candidates.size() == 3'975'589);

    broad_phase->clear();
}

#ifdef IPC_TOOLKIT_WITH_CUDA
TEST_CASE("Puffer-Ball", "[ccd][broad_phase][stq][cuda]")
{
    Eigen::MatrixXd V0, V1;
    Eigen::MatrixXi E, F;

    if (!tests::load_mesh("private/puffer-ball/20.ply", V0, E, F)
        || !tests::load_mesh("private/puffer-ball/21.ply", V1, E, F)) {
        SKIP("Puffer-ball meshes not available");
    }

    CollisionMesh mesh(V0, E, F);

    const auto stq = std::make_shared<SweepAndTiniestQueue>();

    Candidates candidates;
    candidates.build(mesh, V0, V1, /*inflation_radius=*/0, stq);

    CHECK(candidates.size() == 249'805'425);
    CHECK(candidates.vv_candidates.size() == 0);
    CHECK(candidates.ev_candidates.size() == 0);
    CHECK(candidates.ee_candidates.size() == 178'227'707);
    CHECK(candidates.fv_candidates.size() == 71'577'718);
}
#endif