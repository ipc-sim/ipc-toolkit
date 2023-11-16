#include <tests/utils.hpp>

#include <catch2/catch_test_macros.hpp>

#include <ipc/ipc.hpp>

using namespace ipc;

TEST_CASE("Compute CFL stepsize", "[ccd][cfl]")
{
    const double dhat = GENERATE(1e-6, 1e-3, 1e-1);

    Eigen::MatrixXd V0, V1;
    Eigen::MatrixXi E, F;
    SECTION("Cube-Cube")
    {
        Eigen::MatrixXi E1, F1;
        const bool success = tests::load_mesh("two-cubes-close.obj", V0, E, F)
            && tests::load_mesh("two-cubes-intersecting.obj", V1, E1, F1);
        REQUIRE(success);
    }
#ifdef NDEBUG
    SECTION("Cloth-Ball")
    {
        const bool success = tests::load_mesh("cloth_ball92.ply", V0, E, F)
            && tests::load_mesh("cloth_ball93.ply", V1, E, F);
        REQUIRE(success);
    }
#endif

    CollisionMesh mesh(V0, E, F);

    Candidates candidates;
    candidates.build(mesh, V0, dhat / 2);

    const double cfl_step_size =
        candidates.compute_cfl_stepsize(mesh, V0, V1, dhat);

    const double ccd_step_size =
        ipc::compute_collision_free_stepsize(mesh, V0, V1);

    CHECK(cfl_step_size <= ccd_step_size);
}
