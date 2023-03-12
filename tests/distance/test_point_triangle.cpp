#include <catch2/catch_all.hpp>

#include <finitediff.hpp>
#include <igl/PI.h>

#include <ipc/distance/point_triangle.hpp>
#include <ipc/utils/eigen_ext.hpp>

using namespace ipc;

inline Eigen::Vector2d
edge_normal(const Eigen::Vector2d& e0, const Eigen::Vector2d& e1)
{
    Eigen::Vector2d e = e1 - e0;
    Eigen::Vector2d normal(-e.y(), e.x());
    return normal.normalized();
}

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
        Eigen::Vector2d perp = edge_normal(
            Eigen::Vector2d(t0.x(), t0.z()), Eigen::Vector2d(t1.x(), t1.z()));
        double scale = GENERATE(0, 1e-12, 1e-4, 1, 2, 11, 1000);
        p.x() = closest_point.x() + scale * perp.x();
        p.z() = closest_point.z() + scale * perp.y();
    }
    SECTION("closest to t1t2")
    {
        double alpha = GENERATE(0.0, 1e-4, 0.5, 1.0 - 1e-4, 1.0);
        closest_point = (t2 - t1) * alpha + t1;
        Eigen::Vector2d perp = edge_normal(
            Eigen::Vector2d(t1.x(), t1.z()), Eigen::Vector2d(t2.x(), t2.z()));
        double scale = GENERATE(0, 1e-12, 1e-4, 1, 2, 11, 1000);
        p.x() = closest_point.x() + scale * perp.x();
        p.z() = closest_point.z() + scale * perp.y();
    }
    SECTION("closest to t2t0")
    {
        double alpha = GENERATE(0.0, 1e-4, 0.5, 1.0 - 1e-4, 1.0);
        closest_point = (t0 - t2) * alpha + t2;
        Eigen::Vector2d perp = edge_normal(
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
        Eigen::Vector2d perp = edge_normal(
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
        Eigen::Vector2d perp = edge_normal(
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
        Eigen::Vector2d perp = edge_normal(
            Eigen::Vector2d(t2.x(), t2.z()), Eigen::Vector2d(t0.x(), t0.z()));
        // double scale = GENERATE(0, 1e-12, 1e-4, 1, 2, 11, 1000);
        double scale = GENERATE(0, 1e-12, 1e-4, 1, 2, 11);
        p.x() = closest_point.x() + scale * perp.x();
        p.z() = closest_point.z() + scale * perp.y();
    }

    Eigen::VectorXd x(12);
    x.segment<3>(0) = p;
    x.segment<3>(3) = t0;
    x.segment<3>(6) = t1;
    x.segment<3>(9) = t2;
    // Compute the gradient using finite differences
    auto f = [&](const Eigen::VectorXd& x) {
        return point_triangle_distance(
            x.segment<3>(0), x.segment<3>(3), x.segment<3>(6), x.segment<3>(9));
    };

    Eigen::VectorXd grad;
    point_triangle_distance_gradient(p, t0, t1, t2, grad);

    Eigen::VectorXd fgrad;
    fd::finite_gradient(x, f, fgrad);

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
        Eigen::Vector2d perp = edge_normal(
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
        Eigen::Vector2d perp = edge_normal(
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
        Eigen::Vector2d perp = edge_normal(
            Eigen::Vector2d(t2.x(), t2.z()), Eigen::Vector2d(t0.x(), t0.z()));
        // double scale = GENERATE(0, 1e-12, 1e-4, 1, 2, 11, 1000);
        double scale = GENERATE(0, 1e-12, 1e-4, 1, 2, 11);
        p.x() = closest_point.x() + scale * perp.x();
        p.z() = closest_point.z() + scale * perp.y();
    }

    PointTriangleDistanceType dtype =
        point_triangle_distance_type(p, t0, t1, t2);

    Eigen::VectorXd x(12);
    x.segment<3>(0) = p;
    x.segment<3>(3) = t0;
    x.segment<3>(6) = t1;
    x.segment<3>(9) = t2;
    // Compute the gradient using finite differences
    auto f = [&](const Eigen::VectorXd& x) {
        return point_triangle_distance(
            x.segment<3>(0), x.segment<3>(3), x.segment<3>(6), x.segment<3>(9),
            dtype);
    };

    Eigen::MatrixXd hess;
    point_triangle_distance_hessian(p, t0, t1, t2, hess);

    Eigen::MatrixXd fhess;
    fd::finite_hessian(x, f, fhess);

    CAPTURE(dtype);
    CHECK(fd::compare_hessian(hess, fhess, 1e-2));
}
