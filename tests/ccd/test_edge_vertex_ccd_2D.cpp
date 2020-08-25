#include <catch2/catch.hpp>

#include <ccd/edge_vertex_ccd_2D.hpp>

using namespace ccd;

TEST_CASE("Vertex-Vertex Impact", "[ccd]")
{
    Eigen::Vector2d v_t0(1.11111, 0.5);
    Eigen::Vector2d v_t1(0.888889, 0.5);
    Eigen::Vector2d e0_t0(1, 0.5);
    Eigen::Vector2d e1_t0(1, 0.75);
    Eigen::Vector2d e0_t1 = e0_t0;
    Eigen::Vector2d e1_t1 = e1_t0;

    double toi, alpha;
    bool is_collision = compute_edge_vertex_time_of_impact(
        e0_t0, e1_t0, v_t0, e0_t1 - e0_t0, e1_t1 - e1_t0, v_t1 - v_t0, //
        toi, alpha);

    CHECK(is_collision);
    if (is_collision) {
        CHECK(toi == Approx(0.5));
        CHECK(alpha == Approx(0));
    }
}
