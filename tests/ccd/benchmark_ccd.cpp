#include <catch2/catch.hpp>

#include <Eigen/Core>

#include <igl/IO>
#include <igl/edges.h>
#include <igl/Timer.h>

#include <ipc/ipc.hpp>

using namespace ipc;

TEST_CASE("Benchmark earliest toi", "[!benchmark][ccd][earliest_toi]")
{
    Eigen::MatrixXd V0, V1;
    Eigen::MatrixXi E, F;
    Eigen::VectorXi group_ids;

    std::string mesh_path =
        std::string(TEST_DATA_DIR) + "slow-broadphase-ccd/0.obj";
    // std::string(TEST_DATA_DIR) + "slow-broadphase-ccd/s0.obj";
    bool success = igl::read_triangle_mesh(mesh_path, V0, F);
    if (!success) {
        return; // Data is private
    }

    mesh_path = std::string(TEST_DATA_DIR) + "slow-broadphase-ccd/1.obj";
    // mesh_path = std::string(TEST_DATA_DIR) + "slow-broadphase-ccd/s1.obj";
    success = igl::read_triangle_mesh(mesh_path, V1, F);
    if (!success) {
        return; // Data is private
    }

    igl::edges(F, E);

    CollisionMesh mesh = CollisionMesh::build_from_full_mesh(V0, E, F);
    // Discard codimensional/internal vertices
    V0 = mesh.vertices(V0);
    V1 = mesh.vertices(V1);

    // BroadPhaseMethod method = GENERATE(
    //     BroadPhaseMethod::BRUTE_FORCE, BroadPhaseMethod::HASH_GRID,
    //     BroadPhaseMethod::SPATIAL_HASH);

    double tolerance = 1e-6;
    int max_iterations = 1e7;

    BENCHMARK("Earliest ToI BF")
    {
        double stepsize = compute_collision_free_stepsize(
            mesh, V0, V1, BroadPhaseMethod::BRUTE_FORCE, tolerance,
            max_iterations);
        // IPC_LOG(critical("stepsize={}", stepsize));
    };
    BENCHMARK("Earliest ToI HG")
    {
        double stepsize = compute_collision_free_stepsize(
            mesh, V0, V1, static_cast<BroadPhaseMethod>(1), tolerance,
            max_iterations);
        // IPC_LOG(critical("stepsize={}", stepsize));
    };
    BENCHMARK("Earliest ToI SH")
    {
        double stepsize = compute_collision_free_stepsize(
            mesh, V0, V1, static_cast<BroadPhaseMethod>(2), tolerance,
            max_iterations);
        // IPC_LOG(critical("stepsize={}", stepsize));
    };
#ifdef IPC_TOOLKIT_WITH_CUDA
    BENCHMARK("Earliest ToI STQ")
    {
        double stepsize = compute_collision_free_stepsize(
            mesh, V0, V1, BroadPhaseMethod::SWEEP_AND_TINIEST_QUEUE, tolerance,
            max_iterations);
        // IPC_LOG(critical("stepsize={}", stepsize));
    };
#endif
}
