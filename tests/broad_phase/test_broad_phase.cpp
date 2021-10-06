#include <catch2/catch.hpp>

#include <Eigen/Core>

#include <igl/IO>
#include <igl/edges.h>
#include <igl/Timer.h>

#include <ipc/ipc.hpp>
#include <ipc/broad_phase/hash_grid.hpp>
#include <ipc/broad_phase/spatial_hash.hpp>
#include <ipc/ccd/ccd.hpp>
#include <ipc/utils/faces_to_edges.hpp>

#include <broad_phase/brute_force_comparison.hpp>

using namespace ipc;

TEST_CASE("Vertex-Vertex Broad Phase", "[ccd][broad_phase]")
{
    Eigen::MatrixXd V_t0(4, 2);
    V_t0.row(0) << 1.11111, 0.5;  // edge 0 vertex 0
    V_t0.row(1) << 1.11111, 0.75; // edge 0 vertex 1
    V_t0.row(2) << 1, 0.5;        // edge 1 vertex 0
    V_t0.row(3) << 1, 0.75;       // edge 1 vertex 1

    Eigen::MatrixXd V_t1 = V_t0;
    V_t1.row(0) << 0.888889, 0.5;  // edge 0 vertex 0
    V_t1.row(1) << 0.888889, 0.75; // edge 0 vertex 1

    Eigen::MatrixXi E(2, 2);
    E.row(0) << 1, 0;
    E.row(1) << 2, 3;

    bool ignore_internal_vertices = GENERATE(false, true);

    BroadPhaseMethod method = GENERATE(
        BroadPhaseMethod::BRUTE_FORCE, BroadPhaseMethod::HASH_GRID,
        BroadPhaseMethod::SPATIAL_HASH);

    double tolerance = 1e-6;
    int max_iterations = 1e7;

    bool is_valid_step = ipc::is_step_collision_free(
        V_t0, V_t1, E, /*F=*/Eigen::MatrixXi(), method, tolerance,
        max_iterations, ignore_internal_vertices);

    CAPTURE(ignore_internal_vertices);
    CHECK(!is_valid_step);
}

TEST_CASE("Entire 2D Mesh", "[ccd][broad_phase]")
{
    Eigen::MatrixXd V_t0;
    igl::readCSV(std::string(TEST_DATA_DIR) + "V_t0.csv", V_t0);

    Eigen::MatrixXd V_t1;
    igl::readCSV(std::string(TEST_DATA_DIR) + "V_t1.csv", V_t1);

    Eigen::MatrixXi E;
    igl::readCSV(std::string(TEST_DATA_DIR) + "E.csv", E);
    E.array() -= 1; // NOTE: Convert from OBJ format to index

    bool ignore_internal_vertices = true; // GENERATE(false, true);

    BroadPhaseMethod method = GENERATE(
        BroadPhaseMethod::BRUTE_FORCE, BroadPhaseMethod::HASH_GRID
        // , BroadPhaseMethod::SPATIAL_HASH
    );

    double tolerance = 1e-6;
    int max_iterations = 1e7;

    bool is_valid_step;
    SECTION("2D")
    {
        is_valid_step = ipc::is_step_collision_free(
            V_t0.leftCols(2), V_t1.leftCols(2), E, /*F=*/Eigen::MatrixXi(),
            method, tolerance, max_iterations, ignore_internal_vertices);
    }
    SECTION("3D")
    {
        is_valid_step = ipc::is_step_collision_free(
            V_t0, V_t1, E, /*F=*/Eigen::MatrixXi(), method, tolerance,
            max_iterations, ignore_internal_vertices);
    }

    CAPTURE(method);
    CHECK(!is_valid_step);
}

TEST_CASE(
    "Test construct_constraint_set() with codimensional points",
    "[construct_constraint_set][broad_phase]")
{
    const double dhat = 1e-3;
    Eigen::MatrixXd V_rest, V;
    igl::readDMAT(
        std::string(TEST_DATA_DIR) + "codim-points/V_rest.dmat", V_rest);
    igl::readDMAT(std::string(TEST_DATA_DIR) + "codim-points/V.dmat", V);
    Eigen::MatrixXi E, F;
    igl::readDMAT(std::string(TEST_DATA_DIR) + "codim-points/E.dmat", E);
    igl::readDMAT(std::string(TEST_DATA_DIR) + "codim-points/F.dmat", F);
    const Eigen::MatrixXi F2E = faces_to_edges(F, E);

    const bool ignore_internal_vertices = false;

#ifdef NDEBUG
    BroadPhaseMethod method = GENERATE(
        BroadPhaseMethod::BRUTE_FORCE, BroadPhaseMethod::HASH_GRID //,
        // BroadPhaseMethod::SPATIAL_HASH
    );
#else
    BroadPhaseMethod method =
        GENERATE(BroadPhaseMethod::HASH_GRID //,
                                             // BroadPhaseMethod::SPATIAL_HASH
        );
#endif

    Constraints constraint_set;
    construct_constraint_set(
        V_rest, V, E, F, dhat, constraint_set, F2E, /*dmin=*/0, method,
        ignore_internal_vertices);

    CHECK(constraint_set.size() != 0);
}

TEST_CASE(
    "Compare HashGrid/SpatialHash against brute force",
    "[hash_grid][spatial_hash][brute_force]")
{
    using namespace ipc;

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

        SECTION("Without group ids") {}
        SECTION("With group ids")
        {
            group_ids.resize(4);
            group_ids << 0, 0, 1, 1;
        }

        F.resize(0, 3);

        U = Eigen::MatrixXd::Zero(V0.rows(), V0.cols());
        U.col(1).head(2).setConstant(2);
        U.col(1).tail(2).setConstant(-2);
    }
    SECTION("Complex")
    {
        // #ifdef NDEBUG
        //         std::string filename =
        //             GENERATE(std::string("cube.obj"),
        //             std::string("bunny.obj"));
        // #else
        std::string filename = "cube.obj";
        // #endif
        std::string mesh_path = std::string(TEST_DATA_DIR) + filename;
        bool success = igl::read_triangle_mesh(mesh_path, V0, F);
        REQUIRE(success);
        igl::edges(F, E);

        U = Eigen::MatrixXd::Zero(V0.rows(), V0.cols());
        U.col(1).setOnes();
    }

    Candidates candidates;

    // int check_hash_grid = GENERATE(0, 1);
    int check_hash_grid = 1;

    double inflation_radius = 1e-2; // GENERATE(0.0, 1e-4, 1e-3, 1e-2, 1e-1);

    auto can_collide = [&group_ids](size_t vi, size_t vj) {
        return group_ids.size() == 0 || group_ids(vi) != group_ids(vj);
    };

    std::function<void(const Eigen::MatrixXd&)> build_candidates;
    if (check_hash_grid) {
        build_candidates = [&](const Eigen::MatrixXd& V1) {
            candidates.clear();

            HashGrid hashgrid;
            hashgrid.resize(V0, V1, E, inflation_radius);
            hashgrid.addVertices(V0, V1, inflation_radius);
            hashgrid.addEdges(V0, V1, E, inflation_radius);
            hashgrid.addFaces(V0, V1, F, inflation_radius);

            hashgrid.getVertexEdgePairs(
                E, candidates.ev_candidates, can_collide);
            hashgrid.getEdgeEdgePairs(E, candidates.ee_candidates, can_collide);
            hashgrid.getFaceVertexPairs(
                F, candidates.fv_candidates, can_collide);
        };
    } else {
        build_candidates = [&](const Eigen::MatrixXd& V1) {
            candidates.clear();

            SpatialHash sh;
            sh.build(V0, V1, E, F, inflation_radius);
            // TODO: use can_collide
            group_ids.resize(0);
            sh.queryMeshForCandidates(
                V0, V1, E, F, candidates,
                /*queryEV=*/true, /*queryEE=*/true, /*queryFV=*/true);
        };
    }

    for (int i = 0; i < 2; i++) {
        Eigen::MatrixXd V1 = V0 + U;

        build_candidates(V1);

        brute_force_comparison(
            V0, V1, E, F, group_ids, candidates, inflation_radius);

        U.setRandom();
        U *= 3;
    }
}
