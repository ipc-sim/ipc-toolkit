#include <catch2/catch.hpp>

#include <Eigen/Core>

#include <igl/IO>
#include <igl/edges.h>
#include <igl/Timer.h>

#include <ipc/ipc.hpp>
#include <ipc/broad_phase/hash_grid.hpp>
#include <ipc/broad_phase/brute_force.hpp>
#include <ipc/ccd/ccd.hpp>

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

    bool ignore_internal_vertices = true; // GENERATE(false, true);

    BroadPhaseMethod method = GENERATE(
        BroadPhaseMethod::BRUTE_FORCE, BroadPhaseMethod::HASH_GRID,
        BroadPhaseMethod::SPATIAL_HASH);

    double tolerance = 1e-6;
    int max_iterations = 1e7;

    bool is_valid_step;
    SECTION("2D")
    {
        is_valid_step = ipc::is_step_collision_free(
            V_t0.leftCols(2), V_t1.rightCols(2), E, /*F=*/Eigen::MatrixXi(),
            method, tolerance, max_iterations, ignore_internal_vertices);
    }
    SECTION("3D")
    {
        is_valid_step = ipc::is_step_collision_free(
            V_t0, V_t1, E, /*F=*/Eigen::MatrixXi(), method, tolerance,
            max_iterations, ignore_internal_vertices);
    }

    CAPTURE(ignore_internal_vertices);
    CHECK(!is_valid_step);
}

TEST_CASE(
    "Test construct_constraint_set() with codimensional points",
    "[construct_constraint_set][broad_phase]")
{
    double dhat = 0.1;
    // double dhat = 0.00173123;
    Eigen::MatrixXd V_rest, V;
    igl::readDMAT(
        std::string(TEST_DATA_DIR) + "codim-points/V_rest.dmat", V_rest);
    igl::readDMAT(std::string(TEST_DATA_DIR) + "codim-points/V.dmat", V);
    Eigen::MatrixXi E, F;
    igl::readDMAT(std::string(TEST_DATA_DIR) + "codim-points/E.dmat", E);
    igl::readDMAT(std::string(TEST_DATA_DIR) + "codim-points/F.dmat", F);

    bool ignore_internal_vertices = false;

#ifdef NDEBUG
    BroadPhaseMethod method = GENERATE(
        BroadPhaseMethod::BRUTE_FORCE, BroadPhaseMethod::HASH_GRID //,
        // BroadPhaseMethod::SPATIAL_HASH
    );
#else
    BroadPhaseMethod method = GENERATE(BroadPhaseMethod::HASH_GRID //,
                                       // BroadPhaseMethod::SPATIAL_HASH
    );
#endif

    Constraints constraint_set;
    construct_constraint_set(
        V_rest, V, E, F, dhat, constraint_set, /*F2E=*/Eigen::MatrixXi(),
        /*dmin=*/0, method, ignore_internal_vertices);

    CHECK(constraint_set.size() != 0);
}
