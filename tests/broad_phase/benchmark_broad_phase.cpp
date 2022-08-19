#include <catch2/catch.hpp>

#include <Eigen/Core>

#include <igl/IO>
#include <igl/edges.h>
#include <igl/Timer.h>

#include <ipc/ipc.hpp>
#include <ipc/broad_phase/hash_grid.hpp>
#include <ipc/broad_phase/spatial_hash.hpp>
#include <ipc/broad_phase/brute_force.hpp>

#include "test_utils.hpp"

using namespace ipc;

TEST_CASE("Benchmark broad phase", "[!benchmark][broad_phase]")
{
    Eigen::MatrixXd V0, V1;
    Eigen::MatrixXi E, F;

#ifdef NDEBUG
    std::string filename =
        GENERATE(std::string("cube.obj"), std::string("bunny.obj"));
#else
    std::string filename = "cube.obj";
#endif
    std::string mesh_path = TEST_DATA_DIR + filename;
    bool success = igl::read_triangle_mesh(mesh_path, V0, F);
    REQUIRE(success);

    VectorMax3d min_V = V0.colwise().minCoeff();
    VectorMax3d max_V = V0.colwise().maxCoeff();
    VectorMax3d center = (max_V + min_V) / 2;
    V0.rowwise() -= center.transpose();

    std::string testcase_name;
    SECTION("Squish")
    {
        double s = GENERATE(0.1, 0.0, -0.1);
        V1 = V0;
        V1.col(0) *= s;

        testcase_name = fmt::format("Squish s={}", s);
    }
    SECTION("Interobject")
    {
        VectorMax3d scale = max_V - min_V;
        for (int j = 0; j < V0.cols(); j++) {
            V0.col(j).array() /= scale(j);
        }

        int n = V0.rows();

        Eigen::MatrixXd V(2 * n, V0.cols());
        V.topRows(n) = V0;
        V.block(0, 0, n, 1).array() += 0.5;
        V.bottomRows(n) = V0;
        V.block(n, 0, n, 1).array() -= 0.5;
        V0 = V;

        double dx = GENERATE(0.0, 0.25, 0.75, 1.0);
        V1 = V0;
        V1.block(0, 0, n, 1).array() -= dx;
        V1.block(n, 0, n, 1).array() += dx;

        Eigen::MatrixXi F_stack(2 * F.rows(), F.cols());
        F_stack.topRows(F.rows()) = F;
        F_stack.bottomRows(F.rows()) = F.array() + n;
        F = F_stack;

        testcase_name = fmt::format("Interobject dx={}", dx);
    }

    igl::edges(F, E);

    double inflation_radius = 1e-2; // GENERATE(take(5, random(0.0, 0.1)));

    CollisionMesh mesh = CollisionMesh::build_from_full_mesh(V0, E, F);
    // Discard codimensional/internal vertices
    V0 = mesh.vertices(V0);
    V1 = mesh.vertices(V1);

    std::vector<std::string> BP_names = { "BF", "HG", "SH", "STQ", "GPU_STQ" };
    for (int i = 0; i < NUM_BROAD_PHASE_METHODS; i++) {
        BroadPhaseMethod method = static_cast<BroadPhaseMethod>(i);
        BENCHMARK(fmt::format("BP {} ({})", testcase_name, BP_names[i]))
        {
            Candidates candidates;
            construct_collision_candidates(
                mesh, V0, V1, candidates, inflation_radius, method);
        };
    }
}

TEST_CASE(
    "Benchmark broad phase on real data",
    "[!benchmark][broad_phase][real_data]")
{
    Eigen::MatrixXd V0, V1;
    Eigen::MatrixXi E, F;

    std::string mesh_path_t0, mesh_path_t1;
    SECTION("Data 0")
    {
        mesh_path_t0 = TEST_DATA_DIR + "slow-broadphase-ccd/0.obj";
        mesh_path_t1 = TEST_DATA_DIR + "slow-broadphase-ccd/1.obj";
    }
    SECTION("Cloth-Ball")
    {
        mesh_path_t0 = TEST_DATA_DIR + "cloth_ball92.ply";
        mesh_path_t1 = TEST_DATA_DIR + "cloth_ball93.ply";
    }
    SECTION("Data 1")
    {
        mesh_path_t0 = TEST_DATA_DIR + "slow-broadphase-ccd/s0.obj";
        mesh_path_t1 = TEST_DATA_DIR + "slow-broadphase-ccd/s1.obj";
    }

    if (!igl::read_triangle_mesh(mesh_path_t0, V0, F)
        || !igl::read_triangle_mesh(mesh_path_t1, V1, F)) {
        return; // Data is private
    }

    igl::edges(F, E);

    double inflation_radius = 1e-2; // GENERATE(take(5, random(0.0, 0.1)));

    CollisionMesh mesh = CollisionMesh::build_from_full_mesh(V0, E, F);
    // Discard codimensional/internal vertices
    V0 = mesh.vertices(V0);
    V1 = mesh.vertices(V1);

    std::vector<std::string> BP_names = { "BF", "HG", "SH", "STQ", "GPU_STQ" };
    for (int i = 0; i < NUM_BROAD_PHASE_METHODS; i++) {
        BroadPhaseMethod method = static_cast<BroadPhaseMethod>(i);
        BENCHMARK(fmt::format("BP Real Data ({})", BP_names[i]))
        {
            Candidates candidates;
            construct_collision_candidates(
                mesh, V0, V1, candidates, inflation_radius, method);
        };
    }
}
