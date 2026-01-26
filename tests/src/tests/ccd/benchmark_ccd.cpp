#include <tests/config.hpp>
#include <tests/utils.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>

#include <ipc/ipc.hpp>
#include <ipc/broad_phase/lbvh.hpp>
#include <ipc/ccd/tight_inclusion_ccd.hpp>

using namespace ipc;

TEST_CASE("Benchmark earliest toi", "[!benchmark][ccd][earliest_toi]")
{
    Eigen::MatrixXd V0, V1;
    Eigen::MatrixXi E, F;

    std::string mesh_name_t0, mesh_name_t1;
    SECTION("Data 0")
    {
        mesh_name_t0 = "private/slow-broadphase-ccd/0.ply";
        mesh_name_t1 = "private/slow-broadphase-ccd/1.ply";
    }
    SECTION("Data 1")
    {
        mesh_name_t0 = "private/slow-broadphase-ccd/s0.ply";
        mesh_name_t1 = "private/slow-broadphase-ccd/s1.ply";
    }
    SECTION("Cloth-ball")
    {
        mesh_name_t0 = "cloth_ball92.ply";
        mesh_name_t1 = "cloth_ball93.ply";
    }
    SECTION("Puffer-Ball")
    {
        mesh_name_t0 = "puffer-ball/20.ply";
        mesh_name_t1 = "puffer-ball/21.ply";
    }

    if (!tests::load_mesh(mesh_name_t0, V0, E, F)
        || !tests::load_mesh(mesh_name_t1, V1, E, F)) {
        SKIP("Slow broadphase CCD meshes are private");
    }

    CollisionMesh mesh = CollisionMesh::build_from_full_mesh(V0, E, F);
    // Discard codimensional/internal vertices
    V0 = mesh.vertices(V0);
    V1 = mesh.vertices(V1);

    TightInclusionCCD ccd(/*tolerance=*/1e-6, /*max_iterations=*/1e7);

    LBVH broad_phase;

    Candidates candidates;
    BENCHMARK("Broad-Phase Candidate Generation")
    {
        candidates.build(mesh, V0, V1, 0, &broad_phase);
    };

    double toi;
    BENCHMARK("Earliest ToI Narrow-Phase")
    {
        toi = candidates.compute_collision_free_stepsize(
            mesh, V0, V1, /*min_distance=*/0, ccd);
    };
}
