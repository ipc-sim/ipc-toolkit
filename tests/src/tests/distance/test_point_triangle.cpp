#include <tests/utils.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <ipc/distance/point_point.hpp>
#include <ipc/distance/point_triangle.hpp>
#include <ipc/utils/eigen_ext.hpp>

#include <finitediff.hpp>

using namespace ipc;

namespace {
double point_triangle_distance_stacked(
    const Eigen::VectorXd& x,
    const PointTriangleDistanceType dtype = PointTriangleDistanceType::AUTO)
{
    assert(x.size() == 12);
    return point_triangle_distance(
        x.segment<3>(0), x.segment<3>(3), x.segment<3>(6), x.segment<3>(9),
        dtype);
}
} // namespace

TEST_CASE("Point-triangle distance", "[distance][point-triangle]")
{
    double py = GENERATE(-10, -1, -1e-12, 0, 1e-12, 1, 10);
    Eigen::Vector3d p(0, py, 0);
    Eigen::Vector3d t0(-1, 0, 1);
    Eigen::Vector3d t1(1, 0, 1);
    Eigen::Vector3d t2(0, 0, -1);

    Eigen::Vector3d closest_point;
    SECTION("closest to triangle")
    {
        double pz = GENERATE(0, -1 + 1e-12, -1, 1, 1 - 1e-12);
        p.z() = pz;
        closest_point = p;
        closest_point.y() = 0;
    }
    SECTION("closest to t0")
    {
        double px = GENERATE(-1, -1 - 1e-12, -11);
        p.x() = px;
        p.z() = t0.z();
        closest_point = t0;
    }
    SECTION("closest to t1")
    {
        double px = GENERATE(1, 1 + 1e-12, 11);
        p.x() = px;
        p.z() = t1.z();
        closest_point = t1;
    }
    SECTION("closest to t2")
    {
        double pz = GENERATE(-1, -1 - 1e-12, -11);
        p.z() = pz;
        closest_point = t2;
    }
    SECTION("closest to t0t1")
    {
        double alpha = GENERATE(0.0, 1e-4, 0.5, 1.0 - 1e-4, 1.0);
        closest_point = (t1 - t0) * alpha + t0;
        Eigen::Vector2d perp = tests::edge_normal(
            Eigen::Vector2d(t0.x(), t0.z()), Eigen::Vector2d(t1.x(), t1.z()));
        double scale = GENERATE(0, 1e-12, 1e-4, 1, 2, 11, 1000);
        p.x() = closest_point.x() + scale * perp.x();
        p.z() = closest_point.z() + scale * perp.y();
    }
    SECTION("closest to t1t2")
    {
        double alpha = GENERATE(0.0, 1e-4, 0.5, 1.0 - 1e-4, 1.0);
        closest_point = (t2 - t1) * alpha + t1;
        Eigen::Vector2d perp = tests::edge_normal(
            Eigen::Vector2d(t1.x(), t1.z()), Eigen::Vector2d(t2.x(), t2.z()));
        double scale = GENERATE(0, 1e-12, 1e-4, 1, 2, 11, 1000);
        p.x() = closest_point.x() + scale * perp.x();
        p.z() = closest_point.z() + scale * perp.y();
    }
    SECTION("closest to t2t0")
    {
        double alpha = GENERATE(0.0, 1e-4, 0.5, 1.0 - 1e-4, 1.0);
        closest_point = (t0 - t2) * alpha + t2;
        Eigen::Vector2d perp = tests::edge_normal(
            Eigen::Vector2d(t2.x(), t2.z()), Eigen::Vector2d(t0.x(), t0.z()));
        double scale = GENERATE(0, 1e-12, 1e-4, 1, 2, 11, 1000);
        p.x() = closest_point.x() + scale * perp.x();
        p.z() = closest_point.z() + scale * perp.y();
    }

    double distance = point_triangle_distance(p, t0, t1, t2);
    CAPTURE(py, closest_point.x(), closest_point.y(), closest_point.z());
    CHECK(
        distance
        == Catch::Approx(point_point_distance(p, closest_point)).margin(1e-12));
}

TEST_CASE(
    "Point-triangle distance gradient", "[distance][point-triangle][grad]")
{
    double py = GENERATE(-10, -1, -1e-12, 0, 1e-12, 1, 10);
    Eigen::Vector3d p(0, py, 0);
    Eigen::Vector3d t0(-1, 0, 1);
    Eigen::Vector3d t1(1, 0, 1);
    Eigen::Vector3d t2(0, 0, -1);

    SECTION("closest to triangle")
    {
        double pz = GENERATE(0, -1 + 1e-12, -1, 1, 1 - 1e-12);
        p.z() = pz;
    }
    SECTION("closest to t0")
    {
        double px = GENERATE(-1, -1 - 1e-12, -11);
        p.x() = px;
        p.z() = t0.z();
    }
    SECTION("closest to t1")
    {
        double px = GENERATE(1, 1 + 1e-12, 11);
        p.x() = px;
        p.z() = t1.z();
    }
    SECTION("closest to t2")
    {
        double pz = GENERATE(-1, -1 - 1e-12, -11);
        p.z() = pz;
    }
    SECTION("closest to t0t1")
    {
        double alpha = GENERATE(0.0, 1e-4, 0.5, 1.0 - 1e-4, 1.0);
        Eigen::Vector3d closest_point = (t1 - t0) * alpha + t0;
        Eigen::Vector2d perp = tests::edge_normal(
            Eigen::Vector2d(t0.x(), t0.z()), Eigen::Vector2d(t1.x(), t1.z()));
        // double scale = GENERATE(0, 1e-12, 1e-4, 1, 2, 11, 1000);
        double scale = GENERATE(0, 1e-12, 1e-4, 1, 2, 11);
        p.x() = closest_point.x() + scale * perp.x();
        p.z() = closest_point.z() + scale * perp.y();
    }
    SECTION("closest to t1t2")
    {
        double alpha = GENERATE(0.0, 1e-4, 0.5, 1.0 - 1e-4, 1.0);
        Eigen::Vector3d closest_point = (t2 - t1) * alpha + t1;
        Eigen::Vector2d perp = tests::edge_normal(
            Eigen::Vector2d(t1.x(), t1.z()), Eigen::Vector2d(t2.x(), t2.z()));
        // double scale = GENERATE(0, 1e-12, 1e-4, 1, 2, 11, 1000);
        double scale = GENERATE(0, 1e-12, 1e-4, 1, 2, 11);
        p.x() = closest_point.x() + scale * perp.x();
        p.z() = closest_point.z() + scale * perp.y();
    }
    SECTION("closest to t2t0")
    {
        double alpha = GENERATE(0.0, 1e-4, 0.5, 1.0 - 1e-4, 1.0);
        Eigen::Vector3d closest_point = (t0 - t2) * alpha + t2;
        Eigen::Vector2d perp = tests::edge_normal(
            Eigen::Vector2d(t2.x(), t2.z()), Eigen::Vector2d(t0.x(), t0.z()));
        // double scale = GENERATE(0, 1e-12, 1e-4, 1, 2, 11, 1000);
        double scale = GENERATE(0, 1e-12, 1e-4, 1, 2, 11);
        p.x() = closest_point.x() + scale * perp.x();
        p.z() = closest_point.z() + scale * perp.y();
    }

    Vector12d x;
    x << p, t0, t1, t2;

    const Vector12d grad = point_triangle_distance_gradient(p, t0, t1, t2);

    // Compute the gradient using finite differences
    Eigen::VectorXd fgrad;
    fd::finite_gradient(
        x,
        [](const Eigen::VectorXd& _x) {
            return point_triangle_distance_stacked(_x);
        },
        fgrad);

    CHECK(fd::compare_gradient(grad, fgrad));
}

TEST_CASE("Point-triangle distance hessian", "[distance][point-triangle][hess]")
{
    double py = GENERATE(-10, -1, -1e-12, 0, 1e-12, 1, 10);
    Eigen::Vector3d p(0, py, 0);
    Eigen::Vector3d t0(-1, 0, 1);
    Eigen::Vector3d t1(1, 0, 1);
    Eigen::Vector3d t2(0, 0, -1);

    SECTION("closest to triangle")
    {
        double pz = GENERATE(0, -1 + 1e-12, -1, 1, 1 - 1e-12);
        p.z() = pz;
    }
    SECTION("closest to t0")
    {
        double px = GENERATE(-1, -1 - 1e-12, -11);
        p.x() = px;
        p.z() = t0.z();
    }
    SECTION("closest to t1")
    {
        double px = GENERATE(1, 1 + 1e-12, 11);
        p.x() = px;
        p.z() = t1.z();
    }
    SECTION("closest to t2")
    {
        double pz = GENERATE(-1, -1 - 1e-12, -11);
        p.z() = pz;
    }
    SECTION("closest to t0t1")
    {
        double alpha = GENERATE(0.0, 1e-4, 0.5, 1.0 - 1e-4, 1.0);
        Eigen::Vector3d closest_point = (t1 - t0) * alpha + t0;
        Eigen::Vector2d perp = tests::edge_normal(
            Eigen::Vector2d(t0.x(), t0.z()), Eigen::Vector2d(t1.x(), t1.z()));
        // double scale = GENERATE(0, 1e-12, 1e-4, 1, 2, 11, 1000);
        double scale = GENERATE(0, 1e-12, 1e-4, 1, 2, 11);
        p.x() = closest_point.x() + scale * perp.x();
        p.z() = closest_point.z() + scale * perp.y();
    }
    SECTION("closest to t1t2")
    {
        double alpha = GENERATE(0.0, 1e-4, 0.5, 1.0 - 1e-4, 1.0);
        Eigen::Vector3d closest_point = (t2 - t1) * alpha + t1;
        Eigen::Vector2d perp = tests::edge_normal(
            Eigen::Vector2d(t1.x(), t1.z()), Eigen::Vector2d(t2.x(), t2.z()));
        // double scale = GENERATE(0, 1e-12, 1e-4, 1, 2, 11, 1000);
        double scale = GENERATE(0, 1e-12, 1e-4, 1, 2, 11);
        p.x() = closest_point.x() + scale * perp.x();
        p.z() = closest_point.z() + scale * perp.y();
    }
    SECTION("closest to t2t0")
    {
        double alpha = GENERATE(0.0, 1e-4, 0.5, 1.0 - 1e-4, 1.0);
        Eigen::Vector3d closest_point = (t0 - t2) * alpha + t2;
        Eigen::Vector2d perp = tests::edge_normal(
            Eigen::Vector2d(t2.x(), t2.z()), Eigen::Vector2d(t0.x(), t0.z()));
        // double scale = GENERATE(0, 1e-12, 1e-4, 1, 2, 11, 1000);
        double scale = GENERATE(0, 1e-12, 1e-4, 1, 2, 11);
        p.x() = closest_point.x() + scale * perp.x();
        p.z() = closest_point.z() + scale * perp.y();
    }

    PointTriangleDistanceType dtype =
        point_triangle_distance_type(p, t0, t1, t2);

    Vector12d x;
    x << p, t0, t1, t2;

    const Matrix12d hess = point_triangle_distance_hessian(p, t0, t1, t2);

    // Compute the gradient using finite differences
    Eigen::MatrixXd fhess;
    fd::finite_hessian(
        x,
        [dtype](const Eigen::VectorXd& _x) {
            return point_triangle_distance_stacked(_x, dtype);
        },
        fhess);

    CAPTURE(dtype);
    CHECK(fd::compare_hessian(hess, fhess, 1e-2));
}
