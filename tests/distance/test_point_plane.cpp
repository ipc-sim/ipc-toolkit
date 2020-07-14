#include <catch2/catch.hpp>

#include <distance/point_plane.hpp>

using namespace ipc;

TEST_CASE("Point-plane distance", "[distance][point-plane]")
{
    double x = GENERATE(take(10, random(-100.0, 100.0)));
    double y = GENERATE(take(10, random(-100.0, 100.0)));
    double z = GENERATE(take(10, random(-100.0, 100.0)));
    Eigen::Vector3d p(x, y, z);

    double y_plane = GENERATE(take(10, random(-100.0, 100.0)));
    Eigen::Vector3d t0(-1, y_plane, 0), t1(1, y_plane, -1), t2(1, y_plane, 0);

    double distance = point_plane_distance(p, t0, t1, t2);
    double expected_distance = abs(y - y_plane);
    CHECK(distance == Approx(expected_distance * expected_distance));
}
