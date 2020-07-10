#include <catch2/catch.hpp>

#include <igl/PI.h>

#include <distance/line_line.hpp>
#include <utils/eigen_ext.hpp>

using namespace ipc;

TEST_CASE("Line-line distance", "[distance][line-line]")
{
    double ya = GENERATE(take(10, random(-100.0, 100.0)));
    Eigen::Vector3d ea0(-1, ya, 0), ea1(1, ya, 0);

    double yb = GENERATE(take(10, random(-100.0, 100.0)));
    Eigen::Vector3d eb0(0, yb, -1), eb1(0, yb, 1);

    double distance = line_line_distance(ea0, ea1, eb0, eb1);
    double expected_distance = abs(ya - yb);
#ifdef USE_DISTANCE_SQUARED
    expected_distance *= expected_distance;
#endif
    CHECK(distance == Approx(expected_distance));
}
