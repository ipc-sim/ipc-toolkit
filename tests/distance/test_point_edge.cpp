#include <catch2/catch.hpp>

#include <test_utils.hpp>

#include <finitediff.hpp>
#include <igl/PI.h>

#include <ipc/distance/point_edge.hpp>
#include <ipc/utils/eigen_ext.hpp>

using namespace ipc;

TEST_CASE("Point-edge distance", "[distance][point-edge]")
{
    int dim = GENERATE(2, 3);
    double expected_distance = GENERATE(-10, -1, -1e-12, 0, 1e-12, 1, 10);
    VectorMax3d p = VectorMax3d::Zero(dim);
    p.y() = expected_distance;
    VectorMax3d e0 = VectorMax3d::Zero(dim);
    e0.x() = -10;
    VectorMax3d e1 = VectorMax3d::Zero(dim);
    e1.x() = 10;

    double distance = point_edge_distance(p, e0, e1);
    CHECK(distance == Approx(expected_distance * expected_distance));
}

TEST_CASE("Point-edge distance all types", "[distance][point-edge]")
{
    const double alpha = GENERATE(range(-1.0, 2.0, 0.1));
    const double d = GENERATE(range(-10.0, 10.0, 1.0));
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
            n = cross(e1 - e0, Eigen::Vector3d::UnitX()).normalized();
        }

        const VectorMax3d p = ((e1 - e0) * alpha + e0) + d * n;

        const double expected_distance =
            alpha < 0 ? (e0 - p).norm() : (alpha > 1 ? (e1 - p).norm() : d);

        CAPTURE(alpha, expected_distance, dim);

        const double distance = point_edge_distance(p, e0, e1);
        CHECK(
            distance
            == Approx(expected_distance * expected_distance).margin(1e-15));
    }
}

TEST_CASE("Point-edge distance gradient", "[distance][point-edge][gradient]")
{
    int dim = GENERATE(2, 3);
    double expected_distance = GENERATE(-10, -1, -1e-12, 0, 1e-12, 1, 10);
    VectorMax3d p = VectorMax3d::Zero(dim);
    p.y() = expected_distance;
    VectorMax3d e0 = VectorMax3d::Zero(dim);
    e0.x() = -10;
    VectorMax3d e1 = VectorMax3d::Zero(dim);
    e1.x() = 10;

    Eigen::VectorXd grad;
    point_edge_distance_gradient(p, e0, e1, grad);

    // Compute the gradient using finite differences
    Eigen::VectorXd x(3 * dim);
    x.segment(0 * dim, dim) = p;
    x.segment(1 * dim, dim) = e0;
    x.segment(2 * dim, dim) = e1;
    auto f = [&dim](const Eigen::VectorXd& x) {
        return point_edge_distance(
            x.head(dim), x.segment(dim, dim), x.tail(dim));
    };
    Eigen::VectorXd fgrad;
    fd::finite_gradient(x, f, fgrad);

    CHECK(fd::compare_gradient(grad, fgrad));
}

TEST_CASE("Point-edge distance hessian", "[distance][point-edge][hessian]")
{
    int dim = GENERATE(2, 3);
    double expected_distance = GENERATE(-10, -1, -1e-12, 0, 1e-12, 1, 10);
    VectorMax3d p = VectorMax3d::Zero(dim);
    p.y() = expected_distance;
    VectorMax3d e0 = VectorMax3d::Zero(dim);
    e0.x() = -10;
    VectorMax3d e1 = VectorMax3d::Zero(dim);
    e1.x() = 10;

    Eigen::MatrixXd hess;
    point_edge_distance_hessian(p, e0, e1, hess);

    // Compute the gradient using finite differences
    Eigen::VectorXd x(3 * dim);
    x.segment(0 * dim, dim) = p;
    x.segment(1 * dim, dim) = e0;
    x.segment(2 * dim, dim) = e1;
    auto f = [&dim](const Eigen::VectorXd& x) {
        return point_edge_distance(
            x.head(dim), x.segment(dim, dim), x.tail(dim));
    };
    Eigen::MatrixXd fhess;
    fd::finite_hessian(x, f, fhess);

    CHECK(fd::compare_hessian(hess, fhess, 1e-2));
}
