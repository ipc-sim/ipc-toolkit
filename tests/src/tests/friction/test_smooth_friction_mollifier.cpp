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

    double f2 = ipc::df1_x_minus_f1_over_x3(x, epsv_times_h);

    CHECK(f2 * x == Catch::Approx(fd_f2[0]).margin(MARGIN).epsilon(EPSILON));
}
TEST_CASE("Pairwise friction transition gradient", "[friction][mollifier][pairwise][transition]")
{
    static constexpr double EPSILON = 2e-2;
    static constexpr double MARGIN = 1e-6;
    const double epsv_times_h = std::pow(10, GENERATE(range(-8, 0, 1)));

    // Test with various static and kinetic friction coefficients
    const double static_mu = GENERATE(0.5, 0.6, 0.7);
    const double kinetic_mu = GENERATE(0.3, 0.4, 0.5);

    double x;
    x = std::pow(10, GENERATE(range(-8, 0, 1)));

    if (x == 1e-8 && epsv_times_h == 1e-8) {
        return;
    }

    CAPTURE(x, epsv_times_h, static_mu, kinetic_mu);

    Eigen::Matrix<double, 1, 1> X;
    X << x;

    // Check gradient for pairwise friction with transition
    Eigen::VectorXd fd_f1_over_x(1);
    fd::finite_gradient(
        X,
        [&](const Eigen::VectorXd& _X) {
            return ipc::f0_SF_pairwise_transition(_X[0], epsv_times_h, static_mu, kinetic_mu);
        },
        fd_f1_over_x);

    double f1_over_x = ipc::f1_SF_over_x_pairwise_transition(x, epsv_times_h, static_mu, kinetic_mu);

    CHECK(
        f1_over_x * x == Catch::Approx(fd_f1_over_x[0]).margin(MARGIN).epsilon(EPSILON));

    // Check Hessian for pairwise friction with transition
    if (x == epsv_times_h) {
        return;
    }

    Eigen::VectorXd fd_f2(1);
    fd::finite_gradient(
        X,
        [&](const Eigen::VectorXd& _X) {
            return ipc::f1_SF_over_x_pairwise_transition(_X[0], epsv_times_h, static_mu, kinetic_mu);
        },
        fd_f2);

    double f2 = ipc::df1_x_minus_f1_over_x3_pairwise_transition(x, epsv_times_h, static_mu, kinetic_mu);

    CHECK(f2 * x == Catch::Approx(fd_f2[0]).margin(MARGIN).epsilon(EPSILON));
}

TEST_CASE("Pairwise friction transition edge cases", "[friction][mollifier][pairwise][edge][transition]")
{
    static constexpr double EPSILON = 2e-2;
    static constexpr double MARGIN = 1e-6;
    const double epsv_times_h = std::pow(10, GENERATE(range(-8, 0, 1)));

    // Test with extreme values of static and kinetic friction
    const double static_mu = GENERATE(1.0, 0.9);  // Edge cases for static friction
    const double kinetic_mu = GENERATE(0.1, 0.2); // Edge cases for kinetic friction

    double x;
    x = std::pow(10, GENERATE(range(-8, 0, 1)));

    if (x == 1e-8 && epsv_times_h == 1e-8) {
        return;
    }

    CAPTURE(x, epsv_times_h, static_mu, kinetic_mu);

    Eigen::Matrix<double, 1, 1> X;
    X << x;

    // Check gradient for extreme friction coefficient values
    Eigen::VectorXd fd_f1_over_x(1);
    fd::finite_gradient(
        X,
        [&](const Eigen::VectorXd& _X) {
            return ipc::f0_SF_pairwise_transition(_X[0], epsv_times_h, static_mu, kinetic_mu);
        },
        fd_f1_over_x);

    double f1_over_x = ipc::f1_SF_over_x_pairwise_transition(x, epsv_times_h, static_mu, kinetic_mu);

    CHECK(
        f1_over_x * x == Catch::Approx(fd_f1_over_x[0]).margin(MARGIN).epsilon(EPSILON));

    // Check Hessian for extreme friction coefficient values
    if (x == epsv_times_h) {
        return;
    }

    Eigen::VectorXd fd_f2(1);
    fd::finite_gradient(
        X,
        [&](const Eigen::VectorXd& _X) {
            return ipc::f1_SF_over_x_pairwise_transition(_X[0], epsv_times_h, static_mu, kinetic_mu);
        },
        fd_f2);

    double f2 = ipc::df1_x_minus_f1_over_x3_pairwise_transition(x, epsv_times_h, static_mu, kinetic_mu);

    CHECK(f2 * x == Catch::Approx(fd_f2[0]).margin(MARGIN).epsilon(EPSILON));
}