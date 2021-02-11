#include <catch2/catch.hpp>

#include <ipc/ipc.hpp>
#include <ipc/ccd/ccd.hpp>

#include <test_utils.hpp>

TEST_CASE("Vertex-Vertex Impact", "[ccd]")
{
    Eigen::Vector2d v_t0(1.11111, 0.5);
    Eigen::Vector2d v_t1(0.888889, 0.5);
    Eigen::Vector2d e0_t0(1, 0.5);
    Eigen::Vector2d e1_t0(1, 0.75);
    Eigen::Vector2d e0_t1 = e0_t0;
    Eigen::Vector2d e1_t1 = e1_t0;

    double toi;
    bool is_collision = ipc::point_edge_ccd_2D(
        v_t0, e0_t0, e1_t0, v_t1, e0_t1, e1_t1, toi,
        /*conservative_rescaling=*/1);

    REQUIRE(is_collision);
    CHECK(toi == Approx(0.5));
}

TEST_CASE("Repeated CCD", "[ccd][thisone]")
{
    Eigen::MatrixXd V0, V1;
    Eigen::MatrixXi E, F;

    bool success = load_mesh("ccd-failure/0.obj", V0, E, F);
    if (!success) {
        return;
    }
    // REQUIRE(success);

    success = load_mesh("ccd-failure/1.obj", V1, E, F);
    if (!success) {
        return;
    }
    // REQUIRE(success);

    bool has_collisions = !ipc::is_step_collision_free(V0, V1, E, F);

    double stepsize = ipc::compute_collision_free_stepsize(V0, V1, E, F);

    fmt::print(
        "has_collisions={} stepsize={:.17g}\n", has_collisions, stepsize);

    bool has_collisions2 =
        !ipc::is_step_collision_free(V0, (V1 - V0) * stepsize + V0, E, F);

    double stepsize2 = ipc::compute_collision_free_stepsize(
        V0, (V1 - V0) * stepsize + V0, E, F);

    fmt::print(
        "has_collisions2={} stepsize2={:.17g}\n", has_collisions2, stepsize2);
}
