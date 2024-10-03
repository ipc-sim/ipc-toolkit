#include <tests/config.hpp>
#include <tests/broad_phase/brute_force_comparison.hpp>
#include <tests/utils.hpp>

#include <ipc/broad_phase/brute_force.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <igl/readCSV.h>
#include <igl/readDMAT.h>

using namespace ipc;

void test_face_face_broad_phase(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V0,
    const std::optional<Eigen::MatrixXd>& V1,
    const BroadPhaseMethod method,
    double inflation_radius)
{
    // Face-face collisions
    if (mesh.num_faces() == 0 || method == BroadPhaseMethod::BRUTE_FORCE) {
        return;
    }

    std::shared_ptr<BroadPhase> broad_phase =
        BroadPhase::make_broad_phase(method);
    broad_phase->can_vertices_collide = mesh.can_collide;
    if (V1.has_value()) {
        broad_phase->build(
            V0, V1.value(), mesh.edges(), mesh.faces(), inflation_radius);
    } else {
        broad_phase->build(V0, mesh.edges(), mesh.faces(), inflation_radius);
    }
    std::vector<FaceFaceCandidate> ff_candidates;
    broad_phase->detect_face_face_candidates(ff_candidates);

    BruteForce bf;
    bf.can_vertices_collide = mesh.can_collide;
    if (V1.has_value()) {
        bf.build(V0, V1.value(), mesh.edges(), mesh.faces(), inflation_radius);
    } else {
        bf.build(V0, mesh.edges(), mesh.faces(), inflation_radius);
    }
    std::vector<FaceFaceCandidate> bf_ff_candidates;
    bf.detect_face_face_candidates(bf_ff_candidates);

    CHECK(ff_candidates.size() > 0);
    CHECK(ff_candidates.size() == bf_ff_candidates.size());
}

void test_broad_phase(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const BroadPhaseMethod method,
    const bool expect_collision = true,
    const std::string& cached_bf_candidates = "")
{
    CAPTURE(method);
    REQUIRE(V0.rows() == mesh.num_vertices());
    REQUIRE(V1.rows() == mesh.num_vertices());

    double inflation_radius = 0;

    Candidates candidates;
    candidates.build(mesh, V0, V1, inflation_radius, method);

    if (expect_collision) {
        CHECK(!candidates.is_step_collision_free(mesh, V0, V1));
    }

    if (method != BroadPhaseMethod::BRUTE_FORCE) {
        brute_force_comparison(
            mesh, V0, V1, candidates, inflation_radius, cached_bf_candidates);
    }

    // Face-face collisions
    test_face_face_broad_phase(mesh, V0, V1, method, 0);
}

Candidates test_broad_phase(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    BroadPhaseMethod method,
    double inflation_radius,
    const std::string& cached_bf_candidates = "")
{
    CAPTURE(method);
    REQUIRE(V.rows() == mesh.num_vertices());

    Candidates candidates;
    candidates.build(mesh, V, inflation_radius, method);

    if (method != BroadPhaseMethod::BRUTE_FORCE) {
        brute_force_comparison(
            mesh, V, V, candidates, inflation_radius, cached_bf_candidates);
    }

    // Face-face collisions
    test_face_face_broad_phase(mesh, V, std::nullopt, method, inflation_radius);

    return candidates;
}

TEST_CASE("Vertex-Vertex Broad Phase", "[ccd][broad_phase][2D]")
{
    Eigen::MatrixXd V0(4, 2);
    V0.row(0) << 1.11111, 0.5;  // edge 0 vertex 0
    V0.row(1) << 1.11111, 0.75; // edge 0 vertex 1
    V0.row(2) << 1, 0.5;        // edge 1 vertex 0
    V0.row(3) << 1, 0.75;       // edge 1 vertex 1

    Eigen::MatrixXd V1 = V0;
    V1.row(0) << 0.888889, 0.5;  // edge 0 vertex 0
    V1.row(1) << 0.888889, 0.75; // edge 0 vertex 1

    Eigen::MatrixXi E(2, 2);
    E.row(0) << 1, 0;
    E.row(1) << 2, 3;

    CollisionMesh mesh(V0, E, /*F=*/Eigen::MatrixXi());

    const BroadPhaseMethod method = GENERATE_BROAD_PHASE_METHODS();

    test_broad_phase(mesh, V0, V1, method);
}

#if defined(NDEBUG) || !(defined(WIN32) || defined(_WIN32) || defined(__WIN32))
TEST_CASE("Broad Phase: 2D Mesh", "[ccd][broad_phase][2D]")
#else
TEST_CASE("Broad Phase: 2D Mesh", "[ccd][broad_phase][2D][.]")
#endif
{
    Eigen::MatrixXd tmp;
    REQUIRE(igl::readCSV((tests::DATA_DIR / "mesh-2D/V_t0.csv").string(), tmp));
    const Eigen::MatrixXd V0_full = tmp.leftCols(2);

    REQUIRE(igl::readCSV((tests::DATA_DIR / "mesh-2D/V_t1.csv").string(), tmp));
    const Eigen::MatrixXd V1_full = tmp.leftCols(2);

    Eigen::MatrixXi E;
    REQUIRE(igl::readCSV((tests::DATA_DIR / "mesh-2D/E.csv").string(), E));
    E.array() -= 1; // NOTE: Convert from OBJ format to index

    CollisionMesh mesh = CollisionMesh::build_from_full_mesh(
        V0_full, E, /*F=*/Eigen::MatrixXi());

    const Eigen::MatrixXd V0 = mesh.vertices(V0_full);
    const Eigen::MatrixXd V1 = mesh.vertices(V1_full);

    const BroadPhaseMethod method = GENERATE_BROAD_PHASE_METHODS();

    test_broad_phase(mesh, V0, V1, method);
}

TEST_CASE(
    "Build collisions with codimensional points", "[broad_phase][collisions]")
{
    const double dhat = 1e-3;
    Eigen::MatrixXd V_rest, V;
    igl::readDMAT(
        (tests::DATA_DIR / "codim-points/V_rest.dmat").string(), V_rest);
    igl::readDMAT((tests::DATA_DIR / "codim-points/V.dmat").string(), V);
    Eigen::MatrixXi E, F;
    igl::readDMAT((tests::DATA_DIR / "codim-points/E.dmat").string(), E);
    igl::readDMAT((tests::DATA_DIR / "codim-points/F.dmat").string(), F);

    CollisionMesh mesh(V_rest, E, F);
    CHECK(mesh.num_codim_vertices() > 0);

    BroadPhaseMethod method = GENERATE_BROAD_PHASE_METHODS();

    CAPTURE(method);

    const Candidates candidates = test_broad_phase(mesh, V, method, dhat);

    Collisions collisions;
    collisions.build(candidates, mesh, V, dhat);
    CHECK(collisions.size() != 0);
}

TEST_CASE("Compare BP against brute force", "[broad_phase]")
{
    using namespace ipc;

    const BroadPhaseMethod method = GENERATE_BROAD_PHASE_METHODS();

    Eigen::MatrixXd V0, U;
    Eigen::MatrixXi E, F;
    Eigen::VectorXi group_ids;

    SECTION("Simple")
    {
        V0.resize(4, 3);
        V0.row(0) << -1, -1, 0;
        V0.row(1) << 1, -1, 0;
        V0.row(2) << 0, 1, 1;
        V0.row(3) << 0, 1, -1;

        E.resize(2, 2);
        E.row(0) << 0, 1;
        E.row(1) << 2, 3;

        // SECTION("Without group ids") { }
        // SECTION("With group ids")
        // {
        //     group_ids.resize(4);
        //     group_ids << 0, 0, 1, 1;
        // }

        F.resize(0, 3);

        U = Eigen::MatrixXd::Zero(V0.rows(), V0.cols());
        U.col(1).head(2).setConstant(2);
        U.col(1).tail(2).setConstant(-2);
    }
    SECTION("Complex")
    {
        REQUIRE(tests::load_mesh("cube.obj", V0, E, F));
        U = Eigen::MatrixXd::Zero(V0.rows(), V0.cols());
        U.col(1).setOnes();
    }

    CollisionMesh mesh(V0, E, F);
    mesh.can_collide = [&group_ids](size_t vi, size_t vj) {
        return group_ids.size() == 0 || group_ids(vi) != group_ids(vj);
    };

    // double inflation_radius = 1e-2;

    for (int i = 0; i < 2; i++) {
        Eigen::MatrixXd V1 = V0 + U;

        test_broad_phase(mesh, V0, V1, method, false);

        U.setRandom();
        U *= 3;
    }
}

TEST_CASE("Cloth-Ball", "[ccd][broad_phase][cloth-ball][.]")
{
    Eigen::MatrixXd V0, V1;
    Eigen::MatrixXi E, F;

    REQUIRE(tests::load_mesh("cloth_ball92.ply", V0, E, F));
    REQUIRE(tests::load_mesh("cloth_ball93.ply", V1, E, F));

    CollisionMesh mesh(V0, E, F);

    const BroadPhaseMethod method = GENERATE_BROAD_PHASE_METHODS();

    test_broad_phase(
        mesh, V0, V1, method, true,
        (tests::DATA_DIR / "cloth_ball_bf_ccd_candidates.json").string());
}
