#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <ipc/geometry/intersection.hpp>

using namespace ipc;

TEST_CASE("edge_triangle_intersection barycentric", "[intersections]")
{
    const Eigen::Vector3d t0(-1, -1, 0), t1(1, -1, 0), t2(0, 1, 0);

    SECTION("vertical edge through the interior")
    {
        double a, b, t;
        const bool hit = edge_triangle_intersection(
            Eigen::Vector3d(0, -0.2, -1), Eigen::Vector3d(0, -0.2, 1), t0, t1,
            t2, a, b, t);
        CHECK(hit);
        CHECK_THAT(t, Catch::Matchers::WithinAbs(0.5, 1e-9));
        CHECK(a >= 0.0);
        CHECK(b >= 0.0);
        CHECK(a + b <= 1.0 + 1e-9);
    }

    SECTION("edge entirely above the plane misses")
    {
        double a, b, t;
        CHECK(!edge_triangle_intersection(
            Eigen::Vector3d(0, 0, 0.5), Eigen::Vector3d(0, 0, 1.5), t0, t1, t2,
            a, b, t));
    }
}
