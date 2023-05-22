#include <catch2/catch_all.hpp>

#include <Eigen/Core>

#include <igl/IO>
#include <igl/edges.h>
#include <igl/Timer.h>

#include <ipc/ipc.hpp>

#include "test_utils.hpp"

using namespace ipc;

static const std::vector<std::string> BP_names = {
    "BF",
    "HG",
    "SH",
    "STQ",
    "GPU_STQ",
    "BROADMARK_GPU_LBVH",
    "BROADMARK_KD",
    "BROADMARK_GRID",
    "BROADMARK_GRID_PARALLEL",
    "BROADMARK_SAP",
    "BROADMARK_SAP_PARALLEL",
    "BROADMARK_DBVT_D",
    "BROADMARK_DBVT_F",
    "BROADMARK_ISAP",
    "BROADMARK_TRACY",
    "BROADMARK_TRACY_PARALLEL",
    "BROADMARK_GRID_SAP",
    "BROADMARK_CGAL",
    "BROADMARK_GPU_GRID",
    "BROADMARK_GPU_SAP",
};

TEST_CASE("Benchmark earliest toi", "[!benchmark][ccd][earliest_toi]")
{
    Eigen::MatrixXd V0, V1;
    Eigen::MatrixXi E, F;

    std::string mesh_path_t0, mesh_path_t1;
    const std::string dir = TEST_DATA_DIR;
    SECTION("Data 0")
    {
        mesh_path_t0 = dir + "slow-broadphase-ccd/0.obj";
        mesh_path_t1 = dir + "slow-broadphase-ccd/1.obj";
    }
    SECTION("Data 1")
    {
        mesh_path_t0 = dir + "slow-broadphase-ccd/s0.obj";
        mesh_path_t1 = dir + "slow-broadphase-ccd/s1.obj";
    }
    SECTION("Cloth-ball")
    {
        mesh_path_t0 = dir + "cloth_ball92.ply";
        mesh_path_t1 = dir + "cloth_ball93.ply";
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

    double tolerance = 1e-6;
    int max_iterations = 1e7;

    double stepsize;
    int i = static_cast<int>(BroadPhaseMethod::HASH_GRID);
    // for (int i = 0; i < static_cast<int>(BroadPhaseMethod::NUM_METHODS); i++)
    // {
    BroadPhaseMethod method = static_cast<BroadPhaseMethod>(i);

    // Broad phase
    // BENCHMARK(fmt::format("Earliest ToI Broad-Phase {}", BP_names[i]))
    Candidates candidates;
    candidates.build(mesh, V0, V1, /*inflation_radius=*/0, method);
    // };

    BENCHMARK(fmt::format("Earliest ToI Narrow-Phase"))
    {
        stepsize = candidates.compute_collision_free_stepsize(
            mesh, V0, V1, tolerance, max_iterations);
    };
    // }
}
