#include <catch2/catch.hpp>

#include <Eigen/Core>

#include <igl/IO>
#include <igl/edges.h>

#include <ipc/ipc.hpp>
#include <ipc/broad_phase/broad_phase.hpp>
#include <ipc/ccd/ccd.hpp>

#include "brute_force_comparison.hpp"
#include "test_utils.hpp"

using namespace ipc;

void test_broad_phase(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    BroadPhaseMethod method,
    bool expect_collision = true,
    const std::string& cached_bf_candidates = "")
{
    CAPTURE(method);
    REQUIRE(V0.rows() == mesh.num_vertices());
    REQUIRE(V1.rows() == mesh.num_vertices());

    double inflation_radius = 0;

    Candidates candidates;
    construct_collision_candidates(
        mesh, V0, V1, candidates, inflation_radius, method);

    if (expect_collision) {
        CHECK(!ipc::is_step_collision_free(candidates, mesh, V0, V1));
    }

    if (method != BroadPhaseMethod::BRUTE_FORCE) {
        brute_force_comparison(
            mesh, V0, V1, candidates, inflation_radius, cached_bf_candidates);
    }
}

void test_broad_phase(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    BroadPhaseMethod method,
    double inflation_radius,
    const std::string& cached_bf_candidates = "")
{
    CAPTURE(method);
    REQUIRE(V.rows() == mesh.num_vertices());

    Candidates candidates;
    construct_collision_candidates(
        mesh, V, candidates, inflation_radius, method);

    if (method != BroadPhaseMethod::BRUTE_FORCE) {
        brute_force_comparison(
            mesh, V, V, candidates, inflation_radius, cached_bf_candidates);
    }
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

    BroadPhaseMethod method = GENERATE(
        BroadPhaseMethod::BRUTE_FORCE, BroadPhaseMethod::HASH_GRID,
        BroadPhaseMethod::SPATIAL_HASH);

    test_broad_phase(mesh, V0, V1, method);
}

#if defined(NDEBUG) || !(defined(WIN32) || defined(_WIN32) || defined(__WIN32))
TEST_CASE("Entire 2D Mesh", "[ccd][broad_phase][2D]")
#else
TEST_CASE("Entire 2D Mesh", "[ccd][broad_phase][2D][!hide]")
#endif
{
    Eigen::MatrixXd tmp;
    REQUIRE(igl::readCSV(TEST_DATA_DIR + "mesh-2D/V_t0.csv", tmp));
    const Eigen::MatrixXd V0_full = tmp.leftCols(2);

    REQUIRE(igl::readCSV(TEST_DATA_DIR + "mesh-2D/V_t1.csv", tmp));
    const Eigen::MatrixXd V1_full = tmp.leftCols(2);

    Eigen::MatrixXi E;
    REQUIRE(igl::readCSV(TEST_DATA_DIR + "mesh-2D/E.csv", E));
    E.array() -= 1; // NOTE: Convert from OBJ format to index

    CollisionMesh mesh = CollisionMesh::build_from_full_mesh(
        V0_full, E, /*F=*/Eigen::MatrixXi());

    const Eigen::MatrixXd V0 = mesh.vertices(V0_full);
    const Eigen::MatrixXd V1 = mesh.vertices(V1_full);

    BroadPhaseMethod method = GENERATE(
        BroadPhaseMethod::BRUTE_FORCE, BroadPhaseMethod::HASH_GRID,
        BroadPhaseMethod::SPATIAL_HASH);

    test_broad_phase(mesh, V0, V1, method);
}

TEST_CASE(
    "Test constraint_set.build() with codimensional points",
    "[broad_phase][constraints]")
{
    const double dhat = 1e-3;
    Eigen::MatrixXd V_rest, V;
    igl::readDMAT(TEST_DATA_DIR + "codim-points/V_rest.dmat", V_rest);
    igl::readDMAT(TEST_DATA_DIR + "codim-points/V.dmat", V);
    Eigen::MatrixXi E, F;
    igl::readDMAT(TEST_DATA_DIR + "codim-points/E.dmat", E);
    igl::readDMAT(TEST_DATA_DIR + "codim-points/F.dmat", F);

    CollisionMesh mesh(V_rest, E, F);

    BroadPhaseMethod method = GENERATE_BROAD_PHASE_METHODS();

    CAPTURE(method);

    test_broad_phase(mesh, V, method, dhat);

    Constraints constraint_set;
    constraint_set.build(mesh, V, dhat, /*dmin=*/0, method);
    CHECK(constraint_set.size() != 0);
}

TEST_CASE("Compare BP against brute force", "[broad_phase]")
{
    using namespace ipc;

    BroadPhaseMethod method = GENERATE_BROAD_PHASE_METHODS();

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
        std::string filename = "cube.obj";
        std::string mesh_path = TEST_DATA_DIR + filename;
        bool success = igl::read_triangle_mesh(mesh_path, V0, F);
        REQUIRE(success);
        igl::edges(F, E);

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

TEST_CASE("Cloth-Ball", "[ccd][broad_phase][cloth-ball][!hide]")
{
    Eigen::MatrixXd V0, V1;
    Eigen::MatrixXi E, F;

    REQUIRE(igl::read_triangle_mesh(TEST_DATA_DIR + "cloth_ball92.ply", V0, F));
    REQUIRE(igl::read_triangle_mesh(TEST_DATA_DIR + "cloth_ball93.ply", V1, F));

    igl::edges(F, E);

    CollisionMesh mesh(V0, E, F);

    double inflation_radius = 0;

    BroadPhaseMethod method = GENERATE_BROAD_PHASE_METHODS();

    test_broad_phase(
        mesh, V0, V1, method, true,
        TEST_DATA_DIR + "cloth_ball_bf_ccd_candidated.json");
}
