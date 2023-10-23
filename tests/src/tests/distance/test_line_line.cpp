#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>
#include <catch2/generators/catch_generators_random.hpp>

#include <ipc/distance/line_line.hpp>
#include <ipc/utils/eigen_ext.hpp>

#include <finitediff.hpp>

using namespace ipc;

namespace {
double line_line_distance_stacked(const Eigen::VectorXd& x)
{
    assert(x.size() == 12);
    return line_line_distance(
        x.head<3>(), x.segment<3>(3), x.segment<3>(6), x.tail<3>());
}
} // namespace

TEST_CASE("Line-line distance", "[distance][line-line]")
{
    double ya = GENERATE(take(10, random(-100.0, 100.0)));
    Eigen::Vector3d ea0(-1, ya, 0), ea1(1, ya, 0);

    double yb = GENERATE(take(10, random(-100.0, 100.0)));
    Eigen::Vector3d eb0(0, yb, -1), eb1(0, yb, 1);

    double distance = line_line_distance(ea0, ea1, eb0, eb1);
    double expected_distance = std::abs(ya - yb);
    CHECK(distance == Catch::Approx(expected_distance * expected_distance));
}

TEST_CASE("Line-line distance gradient", "[distance][line-line][gradient]")
{
    double ya = GENERATE(take(10, random(-10.0, 10.0)));
    Eigen::Vector3d ea0(-1, ya, 0), ea1(1, ya, 0);

    double yb = GENERATE(take(10, random(-10.0, 10.0)));
    Eigen::Vector3d eb0(0, yb, -1), eb1(0, yb, 1);

    const Vector12d grad = line_line_distance_gradient(ea0, ea1, eb0, eb1);

    Vector12d x;
    x << ea0, ea1, eb0, eb1;
    Eigen::VectorXd expected_grad;
    expected_grad.resize(grad.size());
    fd::finite_gradient(x, line_line_distance_stacked, expected_grad);

    CAPTURE(ya, yb, (grad - expected_grad).norm());
    bool is_grad_correct = fd::compare_gradient(grad, expected_grad);
    CHECK(is_grad_correct);
}

TEST_CASE("Line-line distance hessian", "[distance][line-line][hessian]")
{
    double ya = GENERATE(take(10, random(-10.0, 10.0)));
    Eigen::Vector3d ea0(-1, ya, 0), ea1(1, ya, 0);

    double yb = GENERATE(take(10, random(-10.0, 10.0)));
    Eigen::Vector3d eb0(0, yb, -1), eb1(0, yb, 1);

    const Matrix12d hess = line_line_distance_hessian(ea0, ea1, eb0, eb1);

    Vector12d x;
    x << ea0, ea1, eb0, eb1;
    Eigen::MatrixXd expected_hess;
    expected_hess.resize(hess.rows(), hess.cols());
    fd::finite_hessian(x, line_line_distance_stacked, expected_hess);

    CAPTURE(ya, yb, (hess - expected_hess).norm());
    CHECK(fd::compare_hessian(hess, expected_hess, 1e-2));
}
