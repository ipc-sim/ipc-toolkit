#include <catch2/catch.hpp>

#include <Eigen/Core>

#include <igl/IO>
#include <igl/edges.h>
#include <igl/Timer.h>

#include <ipc/ipc.hpp>

#include "test_utils.hpp"

using namespace ipc;

TEST_CASE("Benchmark earliest toi", "[!benchmark][ccd][earliest_toi]")
{
    Eigen::MatrixXd V0, V1;
    Eigen::MatrixXi E, F;

    std::string mesh_path_t0, mesh_path_t1;
    const std::string dir = TEST_DATA_DIR + "slow-broadphase-ccd/";
    SECTION("Data 0")
    {
        mesh_path_t0 = dir + "0.obj";
        mesh_path_t1 = dir + "1.obj";
    }
    SECTION("Data 1")
    {
        mesh_path_t0 = dir + "s0.obj";
        mesh_path_t1 = dir + "s1.obj";
    }

    bool success = igl::read_triangle_mesh(mesh_path_t0, V0, F);
    if (!success) {
        return; // Data is private
    }

    success = igl::read_triangle_mesh(mesh_path_t1, V1, F);
    if (!success) {
        return; // Data is private
    }

    igl::edges(F, E);

    CollisionMesh mesh = CollisionMesh::build_from_full_mesh(V0, E, F);
    // Discard codimensional/internal vertices
    V0 = mesh.vertices(V0);
    V1 = mesh.vertices(V1);

    std::vector<std::string> BP_names = { "BF", "HG", "SH", "STQ", "GPU_STQ" };

    double tolerance = 1e-6;
    int max_iterations = 1e7;

    for (int i = 0; i < static_cast<int>(BroadPhaseMethod::NUM_METHODS); i++) {
        BroadPhaseMethod method = static_cast<BroadPhaseMethod>(i);
        BENCHMARK(fmt::format("Earliest ToI {}", BP_names[i]))
        {
            double stepsize = compute_collision_free_stepsize(
                mesh, V0, V1, method, tolerance, max_iterations);
            // IPC_LOG(critical("stepsize={}", stepsize));
        };
    }
}
