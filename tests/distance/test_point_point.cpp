#include <catch2/catch.hpp>

#include <finitediff.hpp>

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

TEST_CASE("Point-point distance gradient", "[distance][point-point][gradient]")
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

    Eigen::VectorXd grad;
    point_point_distance_gradient(p0, p1, grad);

    // Compute the gradient using finite differences
    Eigen::VectorXd x(2 * dim);
    x.head(dim) = p0;
    x.tail(dim) = p1;
    auto f = [&dim](const Eigen::VectorXd& x) {
        return point_point_distance(x.head(dim), x.tail(dim));
    };
    Eigen::VectorXd fgrad;
    fd::finite_gradient(x, f, fgrad);

    CHECK(fd::compare_gradient(grad, fgrad));
}

TEST_CASE("Point-point distance hessian", "[distance][point-point][hessian]")
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

    Eigen::MatrixXd hess;
    point_point_distance_hessian(p0, p1, hess);

    // Compute the gradient using finite differences
    Eigen::VectorXd x(2 * dim);
    x.head(dim) = p0;
    x.tail(dim) = p1;
    auto f = [&dim](const Eigen::VectorXd& x) {
        Eigen::VectorXd grad;
        point_point_distance_gradient(x.head(dim), x.tail(dim), grad);
        return grad;
    };
    Eigen::MatrixXd fhess;
    fd::finite_jacobian(x, f, fhess);

    CAPTURE((hess - fhess).squaredNorm());
    CHECK(fd::compare_hessian(hess, fhess));
}
