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

    BENCHMARK("Earliest ToI")
    {
        double stpesize = compute_collision_free_stepsize(
            V0, V1, E, F, /*ignore_codimensional_vertices=*/false);
    };
}
