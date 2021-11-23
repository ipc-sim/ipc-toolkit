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
        std::string(TEST_DATA_DIR) + "slow-broadphase-ccd/s0.obj";
    bool success = igl::read_triangle_mesh(mesh_path, V0, F);
    if (!success) {
        return; // Data is private
    }

    mesh_path = std::string(TEST_DATA_DIR) + "slow-broadphase-ccd/s1.obj";
    success = igl::read_triangle_mesh(mesh_path, V1, F);
    if (!success) {
        return; // Data is private
    }

    igl::edges(F, E);

    BroadPhaseMethod method = GENERATE(
        BroadPhaseMethod::BRUTE_FORCE, BroadPhaseMethod::HASH_GRID,
        BroadPhaseMethod::SPATIAL_HASH);

    double tolerance = 1e-6;
    int max_iterations = 1e7;

    BENCHMARK("Earliest ToI")
    {
        double stpesize = compute_collision_free_stepsize(
            V0, V1, /*codim_V=*/Eigen::VectorXi(), E, F, method, tolerance,
            max_iterations);
    };
}
