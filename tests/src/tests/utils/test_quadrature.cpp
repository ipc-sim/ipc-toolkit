#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <ipc/utils/math.hpp>
#include <ipc/utils/quadrature.hpp>

using namespace ipc;

constexpr int N = 15;

TEST_CASE("Line Quadrature", "[quadrature]")
{
    Eigen::VectorXd pts, weights;
    line_quadrature(N, pts, weights);

    const double a = 10.2, b = 3.4;
    Eigen::VectorXd ys = pts.array() * a + b;
    CHECK(ys.size() == N);

    const double quadrature = weights.dot(ys);
    const double analytic = 0.5 * a + b;
    CHECK(quadrature == Catch::Approx(analytic).margin(1e-12));
}

TEST_CASE("Triangle Quadrature", "[quadrature]")
{
    Eigen::Matrix<double, -1, 2> pts;
    Eigen::VectorXd weights;
    triangle_quadrature(N, pts, weights);

    const double a = 10.2, b = 3.4, c = 1.2;
    Eigen::VectorXd ys = pts.col(0).array() * a + pts.col(1).array() * b + c;
    CHECK(ys.size() == (N * (N + 1) / 2));

    const double quadrature = weights.dot(ys);
    const double analytic = a / 6.0 + b / 6.0 + c / 2.0;
    CHECK(quadrature == Catch::Approx(analytic).margin(1e-12));
}
