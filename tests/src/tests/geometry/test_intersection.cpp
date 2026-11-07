#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <ipc/geometry/intersection.hpp>

#include <cmath>

using namespace ipc;

TEST_CASE("edge_triangle_intersection barycentric", "[intersections]")
{
    const Eigen::Vector3d t0(-1, -1, 0), t1(1, -1, 0), t2(0, 1, 0);

    SECTION("vertical edge through the interior")
    {
        double u, v, t;
        const bool hit = edge_triangle_intersection(
            Eigen::Vector3d(0, -0.2, -1), Eigen::Vector3d(0, -0.2, 1), t0, t1,
            t2, u, v, t);
        CHECK(hit);
        CHECK_THAT(t, Catch::Matchers::WithinAbs(0.5, 1e-9));
        CHECK(u >= 0.0);
        CHECK(v >= 0.0);
        CHECK(u + v <= 1.0 + 1e-9);
    }

    SECTION("recovered point matches the barycentric coordinates")
    {
        const Eigen::Vector3d e0(0.1, -0.3, -2), e1(0.1, -0.3, 1);

        double u, v, t;
        REQUIRE(edge_triangle_intersection(e0, e1, t0, t1, t2, u, v, t));

        // The two parameterizations of the hit point must agree.
        const Eigen::Vector3d p_tri = t0 + u * (t1 - t0) + v * (t2 - t0);
        const Eigen::Vector3d p_edge = e0 + t * (e1 - e0);
        CHECK((p_tri - p_edge).norm() < 1e-12);
        CHECK_THAT(p_edge.z(), Catch::Matchers::WithinAbs(0.0, 1e-12));
    }

    SECTION("edge entirely above the plane misses")
    {
        double u, v, t;
        CHECK(!edge_triangle_intersection(
            Eigen::Vector3d(0, 0, 0.5), Eigen::Vector3d(0, 0, 1.5), t0, t1, t2,
            u, v, t));
        // Outputs are NaN, not indeterminate, on a miss.
        CHECK(std::isnan(u));
        CHECK(std::isnan(v));
        CHECK(std::isnan(t));
    }

    SECTION("edge crosses the plane outside the triangle")
    {
        double u, v, t;
        CHECK(!edge_triangle_intersection(
            Eigen::Vector3d(5, 5, -1), Eigen::Vector3d(5, 5, 1), t0, t1, t2, u,
            v, t));
        CHECK(std::isnan(u));
        CHECK(std::isnan(v));
        CHECK(std::isnan(t));
    }

    SECTION("agrees with is_edge_intersecting_triangle")
    {
        const Eigen::Vector3d e0(0, -0.2, -1), e1(0, -0.2, 1);
        double u, v, t;
        CHECK(
            edge_triangle_intersection(e0, e1, t0, t1, t2, u, v, t)
            == is_edge_intersecting_triangle(e0, e1, t0, t1, t2));
    }
}
