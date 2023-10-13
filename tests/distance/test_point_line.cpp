#include <catch2/catch_all.hpp>

#include <utils.hpp>

#include <iostream>

#include <finitediff.hpp>

#include <ipc/distance/point_line.hpp>
#include <ipc/utils/eigen_ext.hpp>

using namespace ipc;

TEST_CASE("Point-line distance", "[distance][point-line]")
{
    int dim = GENERATE(2, 3);

    double y_point = GENERATE(take(10, random(-100.0, 100.0)));
    VectorMax3d p = VectorMax3d::Zero(dim);
    p.y() = y_point;

    double y_line = GENERATE(take(10, random(-100.0, 100.0)));
    VectorMax3d e0 = VectorMax3d::Zero(dim);
    VectorMax3d e1 = VectorMax3d::Zero(dim);
    e0.x() = -1;
    e0.y() = y_line;
    e1.x() = 1;
    e1.y() = y_line;

    double distance = point_line_distance(p, e0, e1);
    double expected_distance = std::abs(y_point - y_line);
    CHECK(distance == Catch::Approx(expected_distance * expected_distance));
}

TEST_CASE("Point-line distance 2", "[distance][point-line]")
{
    const double alpha = GENERATE(range(-1.0, 2.0, 0.1));
    const double expected_distance = GENERATE(range(-10.0, 10.0, 1.0));
    const int dim = GENERATE(2, 3);
    const int n_random_edges = 20;

    for (int i = 0; i < n_random_edges; i++) {
        VectorMax3d e0, e1, n;
        if (dim == 2) {
            e0 = Eigen::Vector2d::Random();
            e1 = Eigen::Vector2d::Random();
            n = edge_normal(e0, e1);
        } else {
            e0 = Eigen::Vector3d::Random();
            e1 = Eigen::Vector3d::Random();
            n = Eigen::Vector3d(e1 - e0).cross(Eigen::Vector3d::UnitX());
            n.normalize();
        }

        const VectorMax3d p = ((e1 - e0) * alpha + e0) + expected_distance * n;

        CAPTURE(alpha, expected_distance, dim);

        const double distance = point_line_distance(p, e0, e1);
        CHECK(
            distance
            == Catch::Approx(expected_distance * expected_distance)
                   .margin(1e-15));
    }
}

TEST_CASE("Point-line distance gradient", "[distance][point-line][gradient]")
{
    // int dim = GENERATE(2, 3);
    int dim = 2;

    double y_point = GENERATE(take(10, random(-10.0, 10.0)));
    VectorMax3d p = VectorMax3d::Zero(dim);
    p.y() = y_point;

    double y_line = GENERATE(take(10, random(-10.0, 10.0)));
    VectorMax3d e0 = VectorMax3d::Zero(dim);
    VectorMax3d e1 = VectorMax3d::Zero(dim);
    e0.x() = -1;
    e0.y() = y_line;
    e1.x() = 1;
    e1.y() = y_line;

    const VectorMax9d grad = point_line_distance_gradient(p, e0, e1);

    // Compute the gradient using finite differences
    VectorMax9d x(3 * dim);
    x << p, e0, e1;
    auto f = [&dim](const Eigen::VectorXd& x) {
        return point_line_distance(
            x.head(dim), x.segment(dim, dim), x.tail(dim));
    };
    Eigen::VectorXd fgrad;
    fd::finite_gradient(x, f, fgrad);

    CAPTURE(dim, y_point, y_line);
    CHECK(fd::compare_gradient(grad, fgrad));
}

TEST_CASE("Point-line distance hessian", "[distance][point-line][hessian]")
{
    int dim = GENERATE(2, 3);

    double y_point = GENERATE(take(10, random(-10.0, 10.0)));
    VectorMax3d p = VectorMax3d::Zero(dim);
    p.y() = y_point;

    double y_line = GENERATE(take(10, random(-10.0, 10.0)));
    VectorMax3d e0 = VectorMax3d::Zero(dim);
    VectorMax3d e1 = VectorMax3d::Zero(dim);
    e0.x() = -1;
    e0.y() = y_line;
    e1.x() = 1;
    e1.y() = y_line;

    const MatrixMax9d hess = point_line_distance_hessian(p, e0, e1);

    // Compute the gradient using finite differences
    VectorMax9d x(3 * dim);
    x << p, e0, e1;
    auto f = [&dim](const Eigen::VectorXd& x) {
        return point_line_distance(
            x.head(dim), x.segment(dim, dim), x.tail(dim));
    };
    Eigen::MatrixXd fhess;
    fd::finite_hessian(x, f, fhess);

    CAPTURE(dim, y_point, y_line);
    CHECK(fd::compare_hessian(hess, fhess, 1e-2));
}

TEST_CASE(
    "Point-line distance hessian case 1",
    "[distance][point-line][hessian][case1]")
{
    Eigen::Vector3d p(-10.8386, 10, -3.91955);
    Eigen::Vector3d e0(0, 0, -1);
    Eigen::Vector3d e1(-1, 0, 1);

    const MatrixMax9d hess = point_line_distance_hessian(p, e0, e1);

    // Compute the gradient using finite differences
    Vector9d x;
    x << p, e0, e1;
    auto f = [](const Eigen::VectorXd& x) {
        return point_line_distance(x.head(3), x.segment(3, 3), x.tail(3));
    };
    Eigen::MatrixXd fhess;
    fd::finite_hessian(x, f, fhess);

    CHECK(fd::compare_hessian(hess, fhess, 1e-2));
}
