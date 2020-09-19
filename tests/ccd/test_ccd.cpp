#include <catch2/catch.hpp>

#include <ipc/ccd/ccd.hpp>

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
