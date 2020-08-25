#include <catch2/catch.hpp>

#include <distance/distance_type.hpp>

using namespace ipc;

inline Eigen::Vector2d
edge_normal(const Eigen::Vector2d& e0, const Eigen::Vector2d& e1)
{
    Eigen::Vector2d e = e1 - e0;
    Eigen::Vector2d normal(-e.y(), e.x());
    return normal.normalized();
}

TEST_CASE("Point-edge distance type", "[distance][distance-type][point-edge]")
{
    // TODO: Implement tests
}

TEST_CASE(
    "Point-triangle distance type", "[distance][distance-type][point-triangle]")
{
    double py = GENERATE(-10, -1, -1e-12, 0, 1e-12, 1, 10);
    Eigen::Vector3d p(0, py, 0);
    Eigen::Vector3d t0(-1, 0, 1);
    Eigen::Vector3d t1(1, 0, 1);
    Eigen::Vector3d t2(0, 0, -1);

    PointTriangleDistanceType expected_dtype;

    SECTION("closest to triangle")
    {
        double pz = GENERATE(0, -1 + 1e-12, -1, 1, 1 - 1e-12);
        p.z() = pz;
        expected_dtype = PointTriangleDistanceType::P_T;
    }
    SECTION("closest to t0")
    {
        double delta = GENERATE(1e-8, 1e-4, 0.1, 11);
        p.x() = t0.x() - delta;
        p.z() = t0.z() + delta;
        expected_dtype = PointTriangleDistanceType::P_T0;
    }
    SECTION("closest to t1")
    {
        double delta = GENERATE(1e-8, 1e-4, 0.1, 11);
        p.x() = t1.x() + delta;
        p.z() = t1.z() + delta;
        expected_dtype = PointTriangleDistanceType::P_T1;
    }
    SECTION("closest to t2")
    {
        double pz = GENERATE(-1 - 1e-12, -1.1, -11);
        p.z() = pz;
        expected_dtype = PointTriangleDistanceType::P_T2;
    }
    SECTION("closest to t0t1")
    {
        double alpha = GENERATE(0.0, 1e-4, 0.5, 1.0 - 1e-4, 1.0);
        Eigen::Vector3d closest_point = (t1 - t0) * alpha + t0;
        Eigen::Vector2d perp = edge_normal(
            Eigen::Vector2d(t0.x(), t0.z()), Eigen::Vector2d(t1.x(), t1.z()));
        double scale = GENERATE(1e-12, 1e-4, 1, 2, 11, 1000);
        p.x() = closest_point.x() + scale * perp.x();
        p.z() = closest_point.z() + scale * perp.y();
        // Remove the ambiguity at the end points
        PointEdgeDistanceType pe_dtype = point_edge_distance_type(p, t0, t1);
        switch (pe_dtype) {
        case PointEdgeDistanceType::P_E0:
            expected_dtype = PointTriangleDistanceType::P_T0;
            break;
        case PointEdgeDistanceType::P_E1:
            expected_dtype = PointTriangleDistanceType::P_T1;
            break;
        case PointEdgeDistanceType::P_E:
            expected_dtype = PointTriangleDistanceType::P_E0;
            break;
        }
    }
    SECTION("closest to t1t2")
    {
        double alpha = GENERATE(0.0, 1e-4, 0.5, 1.0 - 1e-4, 1.0);
        Eigen::Vector3d closest_point = (t2 - t1) * alpha + t1;
        Eigen::Vector2d perp = edge_normal(
            Eigen::Vector2d(t1.x(), t1.z()), Eigen::Vector2d(t2.x(), t2.z()));
        double scale = GENERATE(1e-12, 1e-4, 1, 2, 11, 1000);
        p.x() = closest_point.x() + scale * perp.x();
        p.z() = closest_point.z() + scale * perp.y();
        // Remove the ambiguity at the end points
        PointEdgeDistanceType pe_dtype = point_edge_distance_type(p, t1, t2);
        switch (pe_dtype) {
        case PointEdgeDistanceType::P_E0:
            expected_dtype = PointTriangleDistanceType::P_T1;
            break;
        case PointEdgeDistanceType::P_E1:
            expected_dtype = PointTriangleDistanceType::P_T2;
            break;
        case PointEdgeDistanceType::P_E:
            expected_dtype = PointTriangleDistanceType::P_E1;
            break;
        }
    }
    SECTION("closest to t2t0")
    {
        double alpha = GENERATE(0.0, 1e-4, 0.5, 1.0 - 1e-4, 1.0);
        Eigen::Vector3d closest_point = (t0 - t2) * alpha + t2;
        Eigen::Vector2d perp = edge_normal(
            Eigen::Vector2d(t2.x(), t2.z()), Eigen::Vector2d(t0.x(), t0.z()));
        double scale = GENERATE(1e-12, 1e-4, 1, 2, 11, 1000);
        p.x() = closest_point.x() + scale * perp.x();
        p.z() = closest_point.z() + scale * perp.y();
        // Remove the ambiguity at the end points
        PointEdgeDistanceType pe_dtype = point_edge_distance_type(p, t2, t0);
        switch (pe_dtype) {
        case PointEdgeDistanceType::P_E0:
            expected_dtype = PointTriangleDistanceType::P_T2;
            break;
        case PointEdgeDistanceType::P_E1:
            expected_dtype = PointTriangleDistanceType::P_T0;
            break;
        case PointEdgeDistanceType::P_E:
            expected_dtype = PointTriangleDistanceType::P_E2;
            break;
        }
    }

    PointTriangleDistanceType dtype =
        point_triangle_distance_type(p, t0, t1, t2);
    CHECK(dtype == expected_dtype);
}

TEST_CASE("Edge-edge distance type", "[distance][distance-type][edge-edge]")
{
    // TODO: Implement tests
}
