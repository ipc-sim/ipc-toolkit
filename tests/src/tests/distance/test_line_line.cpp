#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>
#include <catch2/generators/catch_generators_random.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>

#include <ipc/distance/line_line.hpp>
#include <ipc/utils/eigen_ext.hpp>
#include <ipc/smooth_contact/distance/edge_edge.hpp>

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

TEST_CASE("Line-line closest point pair", "[distance][line-line]")
{
    double ya = GENERATE(take(5, random(-100.0, 100.0)));
    double za = GENERATE(take(5, random(-0.5, 0.5)));
    Eigen::Vector3d ea0(-1, ya, za), ea1(1, ya, 0);

    double yb = GENERATE(take(5, random(-100.0, 100.0)));
    double zb = GENERATE(take(5, random(-0.5, 0.5)));
    Eigen::Vector3d eb0(0, yb, - 1 - zb), eb1(0, yb, 1);

    double expected_distance = line_line_distance(ea0, ea1, eb0, eb1);

    Eigen::Matrix<double, 3, 2> closest_points
        = line_line_closest_point_pairs<double>(ea0, ea1, eb0, eb1);
    double distance = (closest_points.col(0) - closest_points.col(1)).squaredNorm();
    
    CHECK(distance == Catch::Approx(expected_distance));
}

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

TEST_CASE("Line-line closest direction hessian", "[distance][line-line][hessian]")
{
    double ya = GENERATE(take(2, random(-10.0, 10.0)));
    Eigen::Vector3d ea0(-1, ya, 0), ea1(1, ya, 0);

    double yb = GENERATE(take(2, random(-10.0, 10.0)));
    Eigen::Vector3d eb0(0, yb, -1), eb1(0, yb, 1);

    using T = ADHessian<12>;
    DiffScalarBase::setVariableCount(12);
    const auto x = slice_positions<T, 4, 3>((Vector12d() << ea0, ea1, eb0, eb1).finished());
    BENCHMARK("AutoDiff Hessian") {
        line_line_closest_point_direction<T>(x.row(0), x.row(1), x.row(2), x.row(3));
    };
    BENCHMARK("Hessian") {
        line_line_closest_point_direction_hessian(ea0, ea1, eb0, eb1);
    };
}

TEST_CASE("Line-line closest point pairs gradient", "[distance][line-line][gradient]")
{
    double ya = GENERATE(take(10, random(-10.0, 10.0)));
    Eigen::Vector3d ea0(-1, ya, 0), ea1(1, ya, 0);

    double yb = GENERATE(take(10, random(-10.0, 10.0)));
    Eigen::Vector3d eb0(0, yb, -1), eb1(0, yb, 1);

    using T = ADGrad<12>;
    DiffScalarBase::setVariableCount(12);
    const auto x = slice_positions<T, 4, 3>((Vector12d() << ea0, ea1, eb0, eb1).finished());
    auto yAD = line_line_closest_point_pairs<T>(x.row(0), x.row(1), x.row(2), x.row(3));
    auto [y, grad] = line_line_closest_point_pairs_gradient(ea0, ea1, eb0, eb1);
    for (int i = 0; i < yAD.size(); i++)
    {
        REQUIRE((yAD(i).getValue() - y(i)) < 1e-8);
        REQUIRE((yAD(i).getGradient() - grad.row(i).transpose()).norm() < 1e-6 * grad.row(i).norm());
    }

    // BENCHMARK("AutoDiff Gradient") {
    //     line_line_closest_point_pairs<T>(x.row(0), x.row(1), x.row(2), x.row(3));
    // };
    // BENCHMARK("Gradient") {
    //     line_line_closest_point_pairs_gradient(ea0, ea1, eb0, eb1);
    // };
}

TEST_CASE("Line-line closest point pairs hessian", "[distance][line-line][hessian]")
{
    double ya = GENERATE(take(10, random(-10.0, 10.0)));
    Eigen::Vector3d ea0(-1, ya, 0.05), ea1(1, ya, 0);

    double yb = GENERATE(take(10, random(-10.0, 10.0)));
    Eigen::Vector3d eb0(0, yb, -1), eb1(0, yb, 1.02);

    using T = ADHessian<12>;
    DiffScalarBase::setVariableCount(12);
    const auto x = slice_positions<T, 4, 3>((Vector12d() << ea0, ea1, eb0, eb1).finished());
    auto yAD = line_line_closest_point_pairs<T>(x.row(0), x.row(1), x.row(2), x.row(3));
    auto [y, grad, hess] = line_line_closest_point_pairs_hessian(ea0, ea1, eb0, eb1);
    for (int i = 0; i < yAD.size(); i++)
    {
        REQUIRE((yAD(i).getValue() - y(i)) < 1e-8);
        REQUIRE((yAD(i).getGradient() - grad.row(i).transpose()).norm() < 1e-8 * grad.row(i).norm());
        REQUIRE((yAD(i).getHessian() - hess[i]).norm() < 1e-8 * hess[i].norm());
    }

    // BENCHMARK("AutoDiff Hessian") {
    //     line_line_closest_point_pairs<T>(x.row(0), x.row(1), x.row(2), x.row(3));
    // };
    // BENCHMARK("Hessian") {
    //     line_line_closest_point_pairs_hessian(ea0, ea1, eb0, eb1);
    // };
}