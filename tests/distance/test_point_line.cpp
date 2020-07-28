#include <catch2/catch.hpp>

#include <finitediff.hpp>

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
    CHECK(distance == Approx(expected_distance * expected_distance));
}

TEST_CASE("Point-line distance gradient", "[distance][point-line][gradient]")
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

    Eigen::VectorXd grad;
    point_line_distance_gradient(p, e0, e1, grad);

    // Compute the gradient using finite differences
    Eigen::VectorXd x(3 * dim);
    x.segment(0 * dim, dim) = p;
    x.segment(1 * dim, dim) = e0;
    x.segment(2 * dim, dim) = e1;
    auto f = [&dim](const Eigen::VectorXd& x) {
        return point_line_distance(
            x.head(dim), x.segment(dim, dim), x.tail(dim));
    };
    Eigen::VectorXd fgrad;
    fd::finite_gradient(x, f, fgrad);

    CHECK(fd::compare_gradient(grad, fgrad));
}

TEST_CASE("Point-line distance hessian", "[distance][point-line][hessian]")
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

    Eigen::MatrixXd hess;
    point_line_distance_hessian(p, e0, e1, hess);

    // Compute the gradient using finite differences
    Eigen::VectorXd x(3 * dim);
    x.segment(0 * dim, dim) = p;
    x.segment(1 * dim, dim) = e0;
    x.segment(2 * dim, dim) = e1;
    auto f = [&dim](const Eigen::VectorXd& x) {
        return point_line_distance(
            x.head(dim), x.segment(dim, dim), x.tail(dim));
    };
    Eigen::MatrixXd fhess;
    fd::finite_hessian(x, f, fhess);

    CHECK(fd::compare_hessian(hess, fhess));
}
