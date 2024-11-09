#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>
#include <catch2/catch_approx.hpp>

#include <finitediff.hpp>

#include <ipc/adhesion/adhesion.hpp>

using namespace ipc;

TEST_CASE("Normal adhesion derivatives", "[adhesion][normal]")
{
    bool use_dist_sqr = GENERATE(false, true);
    double dhat_p = GENERATE_COPY(range(use_dist_sqr ? -2 : -5, 0));
    dhat_p = pow(10, dhat_p);
    double dhat_a = 2 * dhat_p;
    const double Y = GENERATE(1e3, 4e3, 1e4, 1e5, 4e5);
    const double eps_c = GENERATE(0.002, 0.05, 0.08, 0.5, 0.7, 3.4);

    double a2 = Y * eps_c / (2 * (dhat_p - dhat_a));

    const double log_min_d = use_dist_sqr ? -4 : -6;
    double d = GENERATE_COPY(range(log_min_d, log10(2 * dhat_a), 0.5));
    d = pow(10, d);
    REQUIRE(d >= 1e-8); // finite difference step size
    Eigen::Matrix<double, 1, 1> d_vec;
    d_vec << d;

    if (d == Catch::Approx(dhat_p)) {
        return;
    }

    // Check gradient

    if (use_dist_sqr) {
        d_vec *= d;
        d *= d;
        REQUIRE(d >= 1e-8); // finite difference step size
        dhat_p *= dhat_p;
        dhat_a *= dhat_a;
        a2 /= 2 * dhat_p * (dhat_p + dhat_a);
    }

    Eigen::VectorXd fgrad(1);
    fd::finite_gradient(
        d_vec,
        [&](const Eigen::VectorXd& d) {
            REQUIRE(d[0] >= 0);
            return normal_adhesion_potential(d[0], dhat_p, dhat_a, a2);
        },
        fgrad);

    Eigen::VectorXd grad(1);
    grad << normal_adhesion_potential_first_derivative(d, dhat_p, dhat_a, a2);

    CAPTURE(dhat_p, dhat_a, d, fgrad(0), grad(0), use_dist_sqr);
    CHECK(fd::compare_gradient(fgrad, grad));

    // Check hessian

    if (std::abs(d - dhat_p) < 1e-8 || std::abs(d - dhat_a) < 1e-8) {
        return; // Adhesion is only C1
    }

    fd::finite_gradient(
        d_vec,
        [&](const Eigen::VectorXd& d) {
            REQUIRE(d[0] >= 0);
            return normal_adhesion_potential_first_derivative(
                d[0], dhat_p, dhat_a, a2);
        },
        fgrad);

    grad << normal_adhesion_potential_second_derivative(d, dhat_p, dhat_a, a2);

    CAPTURE(dhat_p, dhat_a, d, fgrad(0), grad(0), use_dist_sqr);
    CHECK(fd::compare_gradient(fgrad, grad));
}

TEST_CASE("Tangential adhesion mollifier", "[adhesion][tangential]")
{
    static constexpr double EPSILON = 2e-2;
    static constexpr double MARGIN = 1e-6;
    const double eps_a = std::pow(10, GENERATE(range(-8, 0, 1)));
    double x;
    // SECTION("x=0") { x = 0; }
    // SECTION("x=gen") {
    x = std::pow(10, GENERATE(range(-8, 0, 1)));
    // }

    if (x == 1e-8 && eps_a == 1e-8) {
        return;
    }

    CAPTURE(x, eps_a);

    Eigen::Matrix<double, 1, 1> X;
    X << x;

    // Check gradient

    Eigen::VectorXd fd_f1(1);
    fd::finite_gradient(
        X,
        [&](const Eigen::VectorXd& _X) {
            return tangential_adhesion_f0(_X[0], eps_a);
        },
        fd_f1);

    CHECK(
        tangential_adhesion_f1(x, eps_a)
        == Catch::Approx(fd_f1[0]).margin(MARGIN).epsilon(EPSILON));

    CHECK(
        tangential_adhesion_f1_over_x(x, eps_a) * x
        == Catch::Approx(fd_f1[0]).margin(MARGIN).epsilon(EPSILON));

    // Check hessian
    if (x == eps_a) {
        return;
    }

    Eigen::VectorXd fd_f2(1);
    fd::finite_gradient(
        X,
        [&](const Eigen::VectorXd& _X) {
            return tangential_adhesion_f1(_X[0], eps_a);
        },
        fd_f2);

    CHECK(
        tangential_adhesion_f2(x, eps_a)
        == Catch::Approx(fd_f2[0]).margin(MARGIN).epsilon(EPSILON));

    double f2 = tangential_adhesion_f2_x_minus_f1_over_x3(x, eps_a);
    f2 *= x * x * x;
    f2 += tangential_adhesion_f1(x, eps_a);
    f2 /= x;

    CHECK(f2 == Catch::Approx(fd_f2[0]).margin(MARGIN).epsilon(EPSILON));
}