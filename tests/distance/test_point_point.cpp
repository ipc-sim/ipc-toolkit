#include <catch2/catch_all.hpp>

#include <finitediff.hpp>

#include <ipc/distance/point_point.hpp>
#include <ipc/utils/eigen_ext.hpp>

using namespace ipc;

TEST_CASE("Point-point distance", "[distance][point-point]")
{
    int dim = GENERATE(2, 3);
    VectorMax3d p0 = VectorMax3d::Zero(dim);
    VectorMax3d p1 = VectorMax3d::Zero(dim);
    double expected_distance = GENERATE(-10, -1, -1e-12, 0, 1e-12, 1, 10);
    SECTION("Aligned with X-axis") { p1(0) = expected_distance; }
    SECTION("Diagonal vector")
    {
        p1.setOnes();
        p1.normalize();
        p1 *= expected_distance;
    }
    double distance = point_point_distance(p0, p1);
    CHECK(distance == Catch::Approx(expected_distance * expected_distance));
}

TEST_CASE("Point-point distance gradient", "[distance][point-point][gradient]")
{
    int dim = GENERATE(2, 3);
    VectorMax3d p0 = VectorMax3d::Zero(dim);
    VectorMax3d p1 = VectorMax3d::Zero(dim);
    double expected_distance = GENERATE(-10, -1, -1e-12, 0, 1e-12, 1, 10);
    SECTION("Aligned with X-axis") { p1(0) = expected_distance; }
    SECTION("Diagonal vector")
    {
        p1.setOnes();
        p1.normalize();
        p1 *= expected_distance;
    }

    const VectorMax6d grad = point_point_distance_gradient(p0, p1);

    // Compute the gradient using finite differences
    VectorMax6d x(2 * dim);
    x << p0, p1;
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
    VectorMax3d p0 = VectorMax3d::Zero(dim);
    VectorMax3d p1 = VectorMax3d::Zero(dim);
    double expected_distance = GENERATE(-10, -1, -1e-12, 0, 1e-12, 1, 10);
    SECTION("Aligned with X-axis") { p1(0) = expected_distance; }
    SECTION("Diagonal vector")
    {
        p1.setOnes();
        p1.normalize();
        p1 *= expected_distance;
    }

    const MatrixMax6d hess = point_point_distance_hessian(p0, p1);

    // Compute the gradient using finite differences
    VectorMax6d x(2 * dim);
    x << p0, p1;
    auto f = [&dim](const Eigen::VectorXd& x) {
        return point_point_distance_gradient(x.head(dim), x.tail(dim));
    };
    Eigen::MatrixXd fhess;
    fd::finite_jacobian(x, f, fhess);

    CAPTURE((hess - fhess).squaredNorm());
    CHECK(fd::compare_hessian(hess, fhess));
}
