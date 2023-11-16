#include <tests/utils.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>
#include <catch2/generators/catch_generators_random.hpp>

#include <ipc/distance/distance_type.hpp>
#include <ipc/distance/point_triangle.hpp>

using namespace ipc;

TEST_CASE("Point-edge distance type", "[distance][distance-type][point-edge]")
{
    const double alpha = GENERATE(range(-1.0, 2.0, 0.1));
    const double distance = GENERATE(range(-10.0, 10.0, 1.0));
    const int dim = GENERATE(2, 3);
    const int n_random_edges = 20;

    for (int i = 0; i < n_random_edges; i++) {
        VectorMax3d e0, e1, n;
        if (dim == 2) {
            e0 = Eigen::Vector2d::Random();
            e1 = Eigen::Vector2d::Random();
            n = tests::edge_normal(e0, e1);
        } else {
            e0 = Eigen::Vector3d::Random();
            e1 = Eigen::Vector3d::Random();
            n = Eigen::Vector3d(e1 - e0).cross(Eigen::Vector3d::UnitX());
            n.normalize();
        }

        const VectorMax3d p = ((e1 - e0) * alpha + e0) + distance * n;

        const PointEdgeDistanceType dtype = point_edge_distance_type(p, e0, e1);

        CAPTURE(alpha, distance, dim, dtype);

        if (std::abs(alpha) < 1e-8) {
            CHECK(
                (dtype == PointEdgeDistanceType::P_E0
                 || dtype == PointEdgeDistanceType::P_E));
        } else if (std::abs(alpha - 1) < 1e-8) {
            CHECK(
                (dtype == PointEdgeDistanceType::P_E1
                 || dtype == PointEdgeDistanceType::P_E));
        } else if (alpha < 0) {
            CHECK(dtype == PointEdgeDistanceType::P_E0);
        } else if (alpha > 1) {
            CHECK(dtype == PointEdgeDistanceType::P_E1);
        } else {
            CHECK(dtype == PointEdgeDistanceType::P_E);
        }
    }
}

struct RandomBarycentricCoordGenerator
    : Catch::Generators::IGenerator<Eigen::Vector3d> {
    Eigen::Vector3d bc;

    RandomBarycentricCoordGenerator() { next(); }

    // via GeneratorUntypedBase:
    // Attempts to move the generator to the next element.
    // Returns true if successful (and thus has another element that can be
    // read)
    bool next() override
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
    Eigen::Vector3d const& get() const override { return bc; }
};

Catch::Generators::GeneratorWrapper<Eigen::Vector3d>
random_barycentric_coord_generator()
{
    return Catch::Generators::GeneratorWrapper<Eigen::Vector3d>(
        Catch::Detail::make_unique<RandomBarycentricCoordGenerator>());
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
        double pz = GENERATE(0, -1 + 1e-12, 1 - 1e-12);
        p.z() = pz;
        expected_dtype = PointTriangleDistanceType::P_T;
    }
    SECTION("random closest to triangle")
    {
        Eigen::Vector3d bc =
            GENERATE(take(100, random_barycentric_coord_generator()));
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
        Eigen::Vector2d perp = tests::edge_normal(
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
        Eigen::Vector2d perp = tests::edge_normal(
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
        Eigen::Vector2d perp = tests::edge_normal(
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
        Eigen::Vector2d perp = tests::edge_normal(
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
        Eigen::Vector2d perp = tests::edge_normal(
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
        Eigen::Vector2d perp = tests::edge_normal(
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
    double alpha = GENERATE(range(-1.0, 2.0, 0.1));
    double s = GENERATE(range(-10.0, 10.0, 1.0));
    if (s == 0)
        return;
    const bool swap_ea = GENERATE(false, true);
    const bool swap_eb = GENERATE(false, true);
    const bool swap_edges = GENERATE(false, true);
    const int n_random_edges = 20;

    const auto sign = [](double x) { return x < 0 ? -1 : 1; };

    for (int i = 0; i < n_random_edges; i++) {
        Eigen::Vector3d ea0 = Eigen::Vector3d::Random();
        Eigen::Vector3d ea1 = Eigen::Vector3d::Random();
        Eigen::Vector3d n =
            (ea1 - ea0).cross(Eigen::Vector3d::UnitX()).normalized();

        Eigen::Vector3d eb0 = ((ea1 - ea0) * alpha + ea0) + s * n;
        if (alpha < 0) {
            alpha = 0;
            n = eb0 - ea0;
            s = n.norm();
            n /= s;
        } else if (alpha > 1) {
            alpha = 1;
            n = eb0 - ea1;
            s = n.norm();
            n /= s;
        }
        REQUIRE(s != 0);

        Eigen::Vector3d eb1 =
            ((ea1 - ea0) * alpha + ea0) + (s + sign(s) * 1) * n;

        if (swap_ea) {
            std::swap(ea0, ea1);
        }
        if (swap_eb) {
            std::swap(eb0, eb1);
        }
        if (swap_edges) {
            std::swap(ea0, eb0);
            std::swap(ea1, eb1);
        }

        // EA0_EB0=0 ✅
        // EA0_EB1=1 ✅
        // EA1_EB0=2 ✅
        // EA1_EB1=3 ✅
        // EA_EB0 =4 ✅
        // EA_EB1 =5 ✅
        // EA0_EB =6 ✅
        // EA1_EB =7 ✅
        // EA_EB  =8 ❌
        EdgeEdgeDistanceType ea0_eb0, ea1_eb0, ea_eb0;
        if (!swap_edges) {
            if (swap_ea && swap_eb) {
                ea0_eb0 = EdgeEdgeDistanceType::EA1_EB1;
                ea1_eb0 = EdgeEdgeDistanceType::EA0_EB1;
                ea_eb0 = EdgeEdgeDistanceType::EA_EB1;
            } else if (swap_ea) {
                ea0_eb0 = EdgeEdgeDistanceType::EA1_EB0;
                ea1_eb0 = EdgeEdgeDistanceType::EA0_EB0;
                ea_eb0 = EdgeEdgeDistanceType::EA_EB0;
            } else if (swap_eb) {
                ea0_eb0 = EdgeEdgeDistanceType::EA0_EB1;
                ea1_eb0 = EdgeEdgeDistanceType::EA1_EB1;
                ea_eb0 = EdgeEdgeDistanceType::EA_EB1;
            } else {
                ea0_eb0 = EdgeEdgeDistanceType::EA0_EB0;
                ea1_eb0 = EdgeEdgeDistanceType::EA1_EB0;
                ea_eb0 = EdgeEdgeDistanceType::EA_EB0;
            }
        } else {
            if (swap_ea && swap_eb) {
                ea0_eb0 = EdgeEdgeDistanceType::EA1_EB1;
                ea1_eb0 = EdgeEdgeDistanceType::EA1_EB0;
                ea_eb0 = EdgeEdgeDistanceType::EA1_EB;
            } else if (swap_ea) {
                ea0_eb0 = EdgeEdgeDistanceType::EA0_EB1;
                ea1_eb0 = EdgeEdgeDistanceType::EA0_EB0;
                ea_eb0 = EdgeEdgeDistanceType::EA0_EB;
            } else if (swap_eb) {
                ea0_eb0 = EdgeEdgeDistanceType::EA1_EB0;
                ea1_eb0 = EdgeEdgeDistanceType::EA1_EB1;
                ea_eb0 = EdgeEdgeDistanceType::EA1_EB;
            } else {
                ea0_eb0 = EdgeEdgeDistanceType::EA0_EB0;
                ea1_eb0 = EdgeEdgeDistanceType::EA0_EB1;
                ea_eb0 = EdgeEdgeDistanceType::EA0_EB;
            }
        }

        const EdgeEdgeDistanceType dtype =
            edge_edge_distance_type(ea0, ea1, eb0, eb1);

        CAPTURE(alpha, s, dtype, swap_ea, swap_eb, swap_edges);

        if (alpha == Catch::Approx(0).margin(1e-15)) {
            CHECK((dtype == ea0_eb0 || dtype == ea_eb0));
        } else if (alpha == Catch::Approx(1).margin(1e-15)) {
            CHECK((dtype == ea1_eb0 || dtype == ea_eb0));
        } else if (alpha < 0) {
            CHECK(dtype == ea0_eb0);
        } else if (alpha > 1) {
            CHECK(dtype == ea1_eb0);
        } else {
            CHECK(dtype == ea_eb0);
        }
    }
}

TEST_CASE(
    "Edge-edge EA_EB distance type", "[distance][distance-type][edge-edge]")
{
    double alpha = GENERATE(take(5, random(0.01, 0.99)));
    double beta = GENERATE(take(5, random(0.01, 0.99)));
    double s = GENERATE(range(-5.0, 5.0, 1.0));
    const bool swap_ea = GENERATE(false, true);
    const bool swap_eb = GENERATE(false, true);
    const bool swap_edges = GENERATE(false, true);
    const int n_random_edges = 20;

    for (int i = 0; i < n_random_edges; i++) {
        Eigen::Vector3d ea1 = Eigen::Vector3d::Random();
        Eigen::Vector3d eb1 = Eigen::Vector3d::Random();

        Eigen::Vector3d ea0 = alpha / (alpha - 1) * ea1;
        Eigen::Vector3d eb0 = beta / (beta - 1) * eb1;

        Eigen::Vector3d n = ea1.cross(eb1);
        REQUIRE(n.norm() != 0);
        n.normalize();

        ea0 += -s / 2 * n;
        ea1 += -s / 2 * n;
        eb0 += s / 2 * n;
        eb1 += s / 2 * n;

        if (swap_ea) {
            std::swap(ea0, ea1);
        }
        if (swap_eb) {
            std::swap(eb0, eb1);
        }
        if (swap_edges) {
            std::swap(ea0, eb0);
            std::swap(ea1, eb1);
        }

        const EdgeEdgeDistanceType dtype =
            edge_edge_distance_type(ea0, ea1, eb0, eb1);

        CAPTURE(alpha, beta, s, dtype, swap_ea, swap_eb, swap_edges);

        // EA_EB  =8 ✅
        CHECK(dtype == EdgeEdgeDistanceType::EA_EB);
    }
}

TEST_CASE(
    "Edge-edge parallel distance type",
    "[distance][distance-type][edge-edge][parallel]")
{
    const double alpha = GENERATE(take(10, random(0.01, 0.99)));
    const double beta = GENERATE(take(10, random(1.01, 1.99)));
    const double s = GENERATE(take(10, random(-5.0, 5.0)));
    const int n_random_edges = 20;

    for (int i = 0; i < n_random_edges; i++) {
        const Eigen::Vector3d ea0 = Eigen::Vector3d::Random();
        const Eigen::Vector3d ea1 = Eigen::Vector3d::Random();
        const Eigen::Vector3d ea = ea1 - ea0;

        Eigen::Vector3d n = (ea1 - ea0).cross(Eigen::Vector3d::UnitX());
        REQUIRE(n.norm() != 0);
        n.normalize();

        const Eigen::Vector3d eb0 = ea0 + alpha * ea + s * n;
        const Eigen::Vector3d eb1 = ea0 + beta * ea + s * n;
        const Eigen::Vector3d eb = eb1 - eb0;

        REQUIRE(
            std::abs(ea.normalized().dot(eb.normalized())) == Catch::Approx(1));
        REQUIRE(ea.cross(eb).norm() == Catch::Approx(0).margin(1e-14));

        const EdgeEdgeDistanceType dtype =
            edge_edge_distance_type(ea0, ea1, eb0, eb1);

        CAPTURE(alpha, beta, s, dtype);
        CHECK(
            (dtype == EdgeEdgeDistanceType::EA_EB0
             || dtype == EdgeEdgeDistanceType::EA1_EB));
    }
}
