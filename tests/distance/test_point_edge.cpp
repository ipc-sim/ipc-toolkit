#include <catch2/catch.hpp>

#include <igl/PI.h>

#include <distance/point_edge.hpp>
#include <utils/eigen_ext.hpp>

using namespace ipc;

TEST_CASE("Point-edge distance", "[distance][point-edge]")
{
    int dim = GENERATE(2, 3);
    double expected_distance = GENERATE(-10, -1, -1e-12, 0, 1e-12, 1, 10);
    Eigen::VectorX3d p = Eigen::VectorX3d::Zero(dim);
    p.y() = expected_distance;
    Eigen::VectorX3d e0 = Eigen::VectorX3d::Zero(dim);
    e0.x() = -10;
    Eigen::VectorX3d e1 = Eigen::VectorX3d::Zero(dim);
    e1.x() = 10;

    double distance = point_edge_distance(p, e0, e1);
#ifdef USE_DISTANCE_SQUARED
    CHECK(distance == Approx(expected_distance * expected_distance));
#else
    CHECK(distance == Approx(abs(expected_distance)));
#endif
}
