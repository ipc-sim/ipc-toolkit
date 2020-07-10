#include <catch2/catch.hpp>

#include <distance/point_line.hpp>

#include <utils/eigen_ext.hpp>

using namespace ipc;

TEST_CASE("Point-line distance", "[distance][point-line]")
{
    int dim = GENERATE(2, 3);

    double y_point = GENERATE(take(10, random(-100.0, 100.0)));
    Eigen::VectorX3d p = Eigen::VectorX3d::Zero(dim);
    p.y() = y_point;

    double y_line = GENERATE(take(10, random(-100.0, 100.0)));
    Eigen::VectorX3d e0 = Eigen::VectorX3d::Zero(dim);
    Eigen::VectorX3d e1 = Eigen::VectorX3d::Zero(dim);
    e0.x() = -1;
    e0.y() = y_line;
    e1.x() = 1;
    e1.y() = y_line;

    double distance = point_line_distance(p, e0, e1);
    double expected_distance = abs(y_point - y_line);
#ifdef USE_DISTANCE_SQUARED
    expected_distance *= expected_distance;
#endif
    CHECK(distance == Approx(expected_distance));
}
