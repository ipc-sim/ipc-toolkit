#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/generators/catch_generators_range.hpp>

#include <ipc/friction/smooth_friction_mollifier.hpp>

#include <finitediff.hpp>

TEST_CASE("Smooth friction gradient", "[friction][mollifier]")
{
    static constexpr double EPSILON = 2e-2;
    static constexpr double MARGIN = 1e-6;
    const double epsv_times_h = std::pow(10, GENERATE(range(-8, 0, 1)));
    double x;
    // SECTION("x=0") { x = 0; }
    // SECTION("x=gen") {
    x = std::pow(10, GENERATE(range(-8, 0, 1)));
    // }

    if (x == 1e-8 && epsv_times_h == 1e-8) {
        return;
    }

    CAPTURE(x, epsv_times_h);

    Eigen::Matrix<double, 1, 1> X;
    X << x;

    // Check gradient

    Eigen::VectorXd fd_f1_over_x(1);
    fd::finite_gradient(
        X,
        [&](const Eigen::VectorXd& _X) {
            return ipc::f0_SF(_X[0], epsv_times_h);
        },
        fd_f1_over_x);
    // fd_f1_over_x /= x;

    double f1_over_x = ipc::f1_SF_over_x(x, epsv_times_h);

    CHECK(
        f1_over_x * x
        == Catch::Approx(fd_f1_over_x[0]).margin(MARGIN).epsilon(EPSILON));

    // Check hessian
    if (x == epsv_times_h) {
        return;
    }
    Eigen::VectorXd fd_f2(1);
    fd::finite_gradient(
        X,
        [&](const Eigen::VectorXd& _X) {
            return ipc::f1_SF_over_x(_X[0], epsv_times_h);
        },
        fd_f2);
    // fd_f2 /= x;

    double f2 = ipc::df1_x_minus_f1_over_x3(x, epsv_times_h);

    CHECK(f2 * x == Catch::Approx(fd_f2[0]).margin(MARGIN).epsilon(EPSILON));
}