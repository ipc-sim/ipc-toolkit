#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/generators/catch_generators_range.hpp>

#include <ipc/friction/smooth_friction_mollifier.hpp>

#include <finitediff.hpp>

using namespace ipc;

TEST_CASE("Smooth friction mollifier", "[friction][mollifier]")
{

    static constexpr double EPSILON = 2e-2;
    static constexpr double MARGIN = 1e-6;
    const double eps_v = std::pow(10, GENERATE(range(-8, 0, 1)));
    double x;
    // SECTION("x=0") { x = 0; }
    // SECTION("x=gen") {
    x = std::pow(10, GENERATE(range(-8, 0, 1)));
    // }

    if (x == 1e-8 && eps_v == 1e-8) {
        return;
    }

    CAPTURE(x, eps_v);

    Eigen::Matrix<double, 1, 1> X;
    X << x;

    // Check gradient

    Eigen::VectorXd fd_f1(1);
    fd::finite_gradient(
        X,
        [&](const Eigen::VectorXd& _X) {
            return smooth_friction_f0(_X[0], eps_v);
        },
        fd_f1);

    CHECK(
        smooth_friction_f1(x, eps_v)
        == Catch::Approx(fd_f1[0]).margin(MARGIN).epsilon(EPSILON));

    CHECK(
        smooth_friction_f1_over_x(x, eps_v) * x
        == Catch::Approx(fd_f1[0]).margin(MARGIN).epsilon(EPSILON));

    // Check hessian
    if (x == eps_v) {
        return;
    }

    Eigen::VectorXd fd_f2(1);
    fd::finite_gradient(
        X,
        [&](const Eigen::VectorXd& _X) {
            return smooth_friction_f1(_X[0], eps_v);
        },
        fd_f2);

    CHECK(
        smooth_friction_f2(x, eps_v)
        == Catch::Approx(fd_f2[0]).margin(MARGIN).epsilon(EPSILON));

    double f2 = smooth_friction_f2_x_minus_f1_over_x3(x, eps_v);
    f2 *= x * x * x;
    f2 += smooth_friction_f1(x, eps_v);
    f2 /= x;

    CHECK(f2 == Catch::Approx(fd_f2[0]).margin(MARGIN).epsilon(EPSILON));
}