#include <tests/utils.hpp>

#include <catch2/catch_test_macros.hpp>

#include <ipc/ipc.hpp>

using namespace ipc;

TEST_CASE("Is step collision free", "[is_step_collision_free]")
{
    Eigen::MatrixXd V0, V1;
    Eigen::MatrixXi E, F;
    REQUIRE(tests::load_mesh("two-cubes-close.obj", V0, E, F));
    REQUIRE(tests::load_mesh("two-cubes-intersecting.obj", V1, E, F));

    CollisionMesh mesh(V0, E, F);

    CHECK(is_step_collision_free(mesh, V0, V0));
    CHECK(!is_step_collision_free(mesh, V0, V1));
}