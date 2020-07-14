#include <catch2/catch.hpp>

#include <distance/point_point.hpp>
#include <utils/eigen_ext.hpp>

using namespace ipc;

TEST_CASE("Point-point distance", "[distance][point-point]")
{
    int dim = GENERATE(2, 3);
    Eigen::VectorX3d p0 = Eigen::VectorX3d::Zero(dim);
    Eigen::VectorX3d p1 = Eigen::VectorX3d::Zero(dim);
    double expected_distance = GENERATE(-10, -1, -1e-12, 0, 1e-12, 1, 10);
    SECTION("Aligned with X-axis") { p1(0) = expected_distance; }
    SECTION("Diagonal vector")
    {
        p1.setOnes();
        p1.normalize();
        p1 *= expected_distance;
    }
    double distance = point_point_distance(p0, p1);
    CHECK(distance == Approx(expected_distance * expected_distance));
}
