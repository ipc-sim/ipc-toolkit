#include <catch2/catch.hpp>

#include <ipc/friction/closest_point.hpp>

using namespace ipc;

TEST_CASE(
    "Point-triangle closest point", "[friction][point-triangle][closest_point]")
{
    Eigen::Vector3d t0(-1, 0, 1), t1(1, 0, 1), t2(0, 0, -1);
    Eigen::Vector2d expected_coords(0.5, 0.5);
    Eigen::Vector3d p =
        t0 + expected_coords[0] * (t1 - t0) + expected_coords[1] * (t2 - t0);
    // p = 1 * t0 + u * t1 - u * t0 + v * t2 - v * t0
    //   = (1 - u - v) * t0 + u * t1 + v * t2
    //   =  w * t0 + u * t1 + v * t2

    Eigen::Vector2d barycentric_coords =
        point_triangle_closest_point(p, t0, t1, t2);
    Eigen::Vector3d p_actual = t0 + barycentric_coords[0] * (t1 - t0)
        + barycentric_coords[1] * (t2 - t0);
    CAPTURE(barycentric_coords);
    CHECK((p - p_actual).norm() == Approx(0).margin(1e-12));
}

TEST_CASE("Edge-edge closest point", "[friction][edge-edge][closest_point]")
{
    Eigen::Vector3d ea0(-1, 0, 0), ea1(1, 0, 0), eb0(0, 0, -1), eb1(0, 0, 1);

    Eigen::Vector2d barycentric_coords =
        edge_edge_closest_point(ea0, ea1, eb0, eb1);
    CAPTURE(barycentric_coords);
    CHECK(barycentric_coords[0] == Approx(0.5));
    CHECK(barycentric_coords[1] == Approx(0.5));
}

TEST_CASE("Point-edge closest point", "[friction][point-edge][closest_point]")
{
    Eigen::Vector3d p(0, 1, 0), e0(-1, 0, 0), e1(1, 0, 0);

    double alpha = point_edge_closest_point(p, e0, e1);
    CHECK(alpha == Approx(0.5));
}
