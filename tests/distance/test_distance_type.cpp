#include <catch2/catch.hpp>

#include <ipc/distance/distance_type.hpp>
#include <ipc/distance/point_triangle.hpp>

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

struct RandomBarycentricCoordGenerator
    : Catch::Generators::IGenerator<Eigen::Vector3d> {
    Eigen::Vector3d bc;

    RandomBarycentricCoordGenerator() { next(); }

    // via GeneratorUntypedBase:
    // Attempts to move the generator to the next element.
    // Returns true if successful (and thus has another element that can be
    // read)
    virtual bool next() override
    {
        double margin = 1e-8;
        bc.head<2>() =
            (Eigen::Array2d::Random() + 1.0) / 2.0 * (1.0 - 2.0 * margin)
            + margin;
        if (bc(0) + bc(1) >= 1) {
            bc(0) = 1 - bc(0);
            bc(1) = 1 - bc(1);
        }
        bc(2) = 1 - bc(0) - bc(1);
        return true;
    }

    // Precondition:
    // The generator is either freshly constructed or the last call to next()
    // returned true
    virtual Eigen::Vector3d const& get() const override { return bc; }
};

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
        double pz = GENERATE(0, -1 + 1e-12, 1 - 1e-12);
        p.z() = pz;
        expected_dtype = PointTriangleDistanceType::P_T;
    }
    SECTION("random closest to triangle")
    {
        Eigen::Vector3d bc = GENERATE(take(
            100,
            Catch::Generators::GeneratorWrapper<Eigen::Vector3d>(
                std::unique_ptr<RandomBarycentricCoordGenerator>(
                    new RandomBarycentricCoordGenerator()))));
        Eigen::Vector3d point_in_plane = bc(0) * t0 + bc(1) * t1 + bc(2) * t2;
        p.x() = point_in_plane.x();
        p.z() = point_in_plane.z();
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
        double alpha = GENERATE(1e-4, 0.5, 1.0 - 1e-4);
        Eigen::Vector3d closest_point = (t1 - t0) * alpha + t0;
        Eigen::Vector2d perp = edge_normal(
            Eigen::Vector2d(t0.x(), t0.z()), Eigen::Vector2d(t1.x(), t1.z()));
        double scale = GENERATE(1e-12, 1e-4, 1, 2, 11, 1000);
        p.x() = closest_point.x() + scale * perp.x();
        p.z() = closest_point.z() + scale * perp.y();
        expected_dtype = PointTriangleDistanceType::P_E0;
    }
    SECTION("random closest to t0t1")
    {
        double alpha = GENERATE(take(100, random(1e-8, 1 - 1e-8)));
        Eigen::Vector3d closest_point = (t1 - t0) * alpha + t0;
        Eigen::Vector2d perp = edge_normal(
            Eigen::Vector2d(t0.x(), t0.z()), Eigen::Vector2d(t1.x(), t1.z()));
        double scale = GENERATE(take(100, random(1e-12, 1e4)));
        p.x() = closest_point.x() + scale * perp.x();
        p.z() = closest_point.z() + scale * perp.y();
        expected_dtype = PointTriangleDistanceType::P_E0;
    }
    SECTION("closest to t1t2")
    {
        double alpha = GENERATE(1e-4, 0.5, 1.0 - 1e-4);
        Eigen::Vector3d closest_point = (t2 - t1) * alpha + t1;
        Eigen::Vector2d perp = edge_normal(
            Eigen::Vector2d(t1.x(), t1.z()), Eigen::Vector2d(t2.x(), t2.z()));
        double scale = GENERATE(1e-12, 1e-4, 1, 2, 11, 1000);
        p.x() = closest_point.x() + scale * perp.x();
        p.z() = closest_point.z() + scale * perp.y();
        expected_dtype = PointTriangleDistanceType::P_E1;
    }
    SECTION("random closest to t1t2")
    {
        double alpha = GENERATE(take(100, random(1e-8, 1 - 1e-8)));
        Eigen::Vector3d closest_point = (t2 - t1) * alpha + t1;
        Eigen::Vector2d perp = edge_normal(
            Eigen::Vector2d(t1.x(), t1.z()), Eigen::Vector2d(t2.x(), t2.z()));
        double scale = GENERATE(take(100, random(1e-12, 1e4)));
        p.x() = closest_point.x() + scale * perp.x();
        p.z() = closest_point.z() + scale * perp.y();
        expected_dtype = PointTriangleDistanceType::P_E1;
    }
    SECTION("closest to t2t0")
    {
        double alpha = GENERATE(1e-4, 0.5, 1.0 - 1e-4);
        Eigen::Vector3d closest_point = (t0 - t2) * alpha + t2;
        Eigen::Vector2d perp = edge_normal(
            Eigen::Vector2d(t2.x(), t2.z()), Eigen::Vector2d(t0.x(), t0.z()));
        double scale = GENERATE(1e-12, 1e-4, 1, 2, 11, 1000);
        p.x() = closest_point.x() + scale * perp.x();
        p.z() = closest_point.z() + scale * perp.y();
        expected_dtype = PointTriangleDistanceType::P_E2;
    }
    SECTION("random closest to t2t0")
    {
        double alpha = GENERATE(take(100, random(1e-8, 1 - 1e-8)));
        Eigen::Vector3d closest_point = (t0 - t2) * alpha + t2;
        Eigen::Vector2d perp = edge_normal(
            Eigen::Vector2d(t2.x(), t2.z()), Eigen::Vector2d(t0.x(), t0.z()));
        double scale = GENERATE(take(100, random(1e-12, 1e4)));
        p.x() = closest_point.x() + scale * perp.x();
        p.z() = closest_point.z() + scale * perp.y();
        expected_dtype = PointTriangleDistanceType::P_E2;
    }

    PointTriangleDistanceType dtype =
        point_triangle_distance_type(p, t0, t1, t2);
    CHECK(dtype == expected_dtype);
}

TEST_CASE(
    "Point-triangle distance type GH Issue",
    "[distance][distance-type][point-triangle][debug]")
{
    Eigen::Vector3d p(0.488166, 0.0132623, 0.289055);
    Eigen::Vector3d v0(0.456476, 0.0526442, 0.260834);
    Eigen::Vector3d v1(0.609111, 0.0595969, 0.275928);
    Eigen::Vector3d v2(0.431262, 0.0508414, 0.255831);

    CHECK(
        point_triangle_distance_type(p, v0, v1, v2)
        == PointTriangleDistanceType::P_E0);
}

TEST_CASE("Edge-edge distance type", "[distance][distance-type][edge-edge]")
{
    // TODO: Implement tests
}
