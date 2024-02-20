#include <tests/config.hpp>

#include <tests/utils.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <ipc/ipc.hpp>

using namespace ipc;

#ifdef IPC_TOOLKIT_WITH_CUDA

TEST_CASE("GPU CCD", "[ccd][gpu]")
{
    Eigen::MatrixXd V0, V1;
    Eigen::MatrixXi E, F;

    std::string mesh_name_t0, mesh_name_t1;
    SECTION("Cloth-Ball")
    {
        mesh_name_t0 = "cloth_ball92.ply";
        mesh_name_t1 = "cloth_ball93.ply";
    }
    // SECTION("Squishy-Ball")
    // {
    //     mesh_name_t0 = "private/puffer-ball/20.ply";
    //     mesh_name_t1 = "private/puffer-ball/21.ply";
    // }

    if (!tests::load_mesh(mesh_name_t0, V0, E, F)
        || !tests::load_mesh(mesh_name_t1, V1, E, F)) {
        return; // Data is private
    }

    CollisionMesh mesh = CollisionMesh::build_from_full_mesh(V0, E, F);
    // Discard codimensional/internal vertices
    V0 = mesh.vertices(V0);
    V1 = mesh.vertices(V1);

    const double tolerance = 1e-6;
    const int max_iterations = 1e7;
    const double min_distance = 0;

    const double toi_cpu = compute_collision_free_stepsize(
        mesh, V0, V1, BroadPhaseMethod::SWEEP_AND_PRUNE, min_distance,
        tolerance, max_iterations);

    // Got this value from running the code
    CHECK(toi_cpu == Catch::Approx(4.76837158203125000e-06));

    const double toi_gpu = compute_collision_free_stepsize(
        mesh, V0, V1, BroadPhaseMethod::SWEEP_AND_TINIEST_QUEUE, min_distance,
        tolerance, max_iterations);

    // Got this value from running the code
    CHECK(toi_gpu == Catch::Approx(3.05175781250000017e-6));
}

#endif