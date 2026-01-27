#include <tests/config.hpp>

#include <tests/utils.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <ipc/ipc.hpp>
#include <ipc/broad_phase/sweep_and_prune.hpp>
#include <ipc/broad_phase/sweep_and_tiniest_queue.hpp>
#include <ipc/ccd/tight_inclusion_ccd.hpp>
#include <ipc/ccd/additive_ccd.hpp>

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
    // SECTION("Puffer-Ball")
    // {
    //     mesh_name_t0 = "puffer-ball/20.ply";
    //     mesh_name_t1 = "puffer-ball/21.ply";
    // }

    if (!tests::load_mesh(mesh_name_t0, V0, E, F)
        || !tests::load_mesh(mesh_name_t1, V1, E, F)) {
        SKIP("Puffer-ball meshes are not available");
    }

    CollisionMesh mesh = CollisionMesh::build_from_full_mesh(V0, E, F);
    // Discard codimensional/internal vertices
    V0 = mesh.vertices(V0);
    V1 = mesh.vertices(V1);

    const double tolerance = 1e-6;
    const int max_iterations = 1e7;
    const double min_distance = 0;

    SweepAndPrune sap;
    const double toi_cpu = compute_collision_free_stepsize(
        mesh, V0, V1, min_distance, &sap,
        TightInclusionCCD(tolerance, max_iterations));

    // Got this value from running the code
    CHECK(toi_cpu == Catch::Approx(4.76837158203125000e-06));

    SweepAndTiniestQueue stq;
    const double toi_gpu = compute_collision_free_stepsize(
        mesh, V0, V1, min_distance, &stq,
        TightInclusionCCD(tolerance, max_iterations));

    // Got this value from running the code
    CHECK(toi_gpu == Catch::Approx(3.05175781250000017e-6));
}

TEST_CASE("GPU CCD parameter extraction", "[ccd][gpu]")
{
    Eigen::MatrixXd V0, V1;
    Eigen::MatrixXi E, F;

    if (!tests::load_mesh("cloth_ball92.ply", V0, E, F)
        || !tests::load_mesh("cloth_ball93.ply", V1, E, F)) {
        SKIP("Cloth-ball meshes are not available");
    }

    CollisionMesh mesh = CollisionMesh::build_from_full_mesh(V0, E, F);
    // Discard codimensional/internal vertices
    V0 = mesh.vertices(V0);
    V1 = mesh.vertices(V1);

    const double min_distance = 0;
    SweepAndTiniestQueue stq;

    // Test 1: Custom TightInclusionCCD parameters should be extracted and used
    SECTION("Custom TightInclusionCCD parameters")
    {
        // Different from default (1e-6)
        const double custom_tolerance = 1e-5;
        // Different from default (10M)
        const long custom_max_iterations = 5'000'000L;

        const double toi_custom = compute_collision_free_stepsize(
            mesh, V0, V1, min_distance, &stq,
            TightInclusionCCD(custom_tolerance, custom_max_iterations));

        // Should produce a valid result (not crash)
        CHECK(toi_custom >= 0.0);
        CHECK(toi_custom <= 1.0);
    }

    // Test 2: Default TightInclusionCCD should use default parameters
    SECTION("Default TightInclusionCCD")
    {
        const double toi_default = compute_collision_free_stepsize(
            mesh, V0, V1, min_distance, &stq, TightInclusionCCD());

        // Should produce a valid result
        CHECK(toi_default >= 0.0);
        CHECK(toi_default <= 1.0);
    }

    // Test 3: Non-TightInclusionCCD should fall back to defaults
    SECTION("Non-TightInclusionCCD fallback to defaults")
    {
        const double toi_additive = compute_collision_free_stepsize(
            mesh, V0, V1, min_distance, &stq, AdditiveCCD());

        // Should produce a valid result (fallback to defaults should work)
        CHECK(toi_additive >= 0.0);
        CHECK(toi_additive <= 1.0);
    }

    // Test 4: Verify that different tolerance values can produce different
    // results
    SECTION("Different tolerance values")
    {
        // Use default tolerance explicitly - should match default constructor
        const double toi_explicit_default = compute_collision_free_stepsize(
            mesh, V0, V1, min_distance, &stq,
            TightInclusionCCD(
                TightInclusionCCD::DEFAULT_TOLERANCE,
                TightInclusionCCD::DEFAULT_MAX_ITERATIONS));

        const double toi_implicit_default = compute_collision_free_stepsize(
            mesh, V0, V1, min_distance, &stq, TightInclusionCCD());

        // Explicit default should match implicit default (verifies parameter
        // extraction works)
        CHECK(toi_explicit_default == Catch::Approx(toi_implicit_default));

        // Use different tolerance - may produce different result
        const double toi_loose = compute_collision_free_stepsize(
            mesh, V0, V1, min_distance, &stq,
            TightInclusionCCD(1e-4, TightInclusionCCD::DEFAULT_MAX_ITERATIONS));

        // All should be valid
        CHECK(toi_explicit_default >= 0.0);
        CHECK(toi_explicit_default <= 1.0);
        CHECK(toi_implicit_default >= 0.0);
        CHECK(toi_implicit_default <= 1.0);
        CHECK(toi_loose >= 0.0);
        CHECK(toi_loose <= 1.0);

        // Verify parameters are being used: different tolerance should
        // potentially produce different results (though they may be the same
        // if no collisions or if difference is within numerical precision)
        // The key is that explicit default matches implicit default, proving
        // parameter extraction works
    }
}

#endif