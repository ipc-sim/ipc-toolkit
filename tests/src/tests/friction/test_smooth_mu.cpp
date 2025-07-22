#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/generators/catch_generators_range.hpp>

#include <ipc/friction/smooth_mu.hpp>

#include <finitediff.hpp>

using namespace ipc;

TEST_CASE("Smooth mu", "[friction][mollifier][mu]")
{

    static constexpr double EPSILON = 1e-4;
    static constexpr double MARGIN = 1e-6;
    // NOTE: Use h=1e-10 for finite difference because min eps_v=1e-8
    static constexpr double H = 1e-10;
    const double mu_s = GENERATE(range(0.0, 1.0, 0.1));
    const double mu_k = GENERATE(range(0.0, 1.0, 0.1));
    const double eps_v = std::pow(10, GENERATE(range(-8, 0, 1)));
    const double x = std::pow(10, GENERATE(range(-8, 0, 1)));

    if (x == 1e-8 && eps_v == 1e-8) {
        return;
    }

    CAPTURE(x, eps_v, x / eps_v, mu_s, mu_k);

    Eigen::Matrix<double, 1, 1> X;
    X << x;

    // Check gradient

    Eigen::VectorXd fd_f1(1);
    fd::finite_gradient(
        X,
        [&](const Eigen::VectorXd& _X) {
            return smooth_mu_f0(_X[0], mu_s, mu_k, eps_v);
        },
        fd_f1, fd::AccuracyOrder::SECOND, H);

    CHECK(
        smooth_mu_f1(x, mu_s, mu_k, eps_v)
        == Catch::Approx(fd_f1[0]).margin(MARGIN).epsilon(EPSILON));

    CHECK(
        smooth_mu_f1_over_x(x, mu_s, mu_k, eps_v) * x
        == Catch::Approx(fd_f1[0]).margin(MARGIN).epsilon(EPSILON));

    // Check hessian
    if (x == eps_v) {
        return;
    }

    Eigen::VectorXd fd_mu1(1);
    fd::finite_gradient(
        X,
        [&](const Eigen::VectorXd& _X) {
            return smooth_mu(_X[0], mu_s, mu_k, eps_v);
        },
        fd_mu1, fd::AccuracyOrder::SECOND, H);

    CHECK(
        smooth_mu_derivative(x, mu_s, mu_k, eps_v)
        == Catch::Approx(fd_mu1[0]).margin(MARGIN).epsilon(EPSILON));

    Eigen::VectorXd fd_f2(1);
    fd::finite_gradient(
        X,
        [&](const Eigen::VectorXd& _X) {
            return smooth_mu_f1(_X[0], mu_s, mu_k, eps_v);
        },
        fd_f2, fd::AccuracyOrder::SECOND, H);

    CHECK(
        smooth_mu_f2(x, mu_s, mu_k, eps_v)
        == Catch::Approx(fd_f2[0]).margin(MARGIN).epsilon(EPSILON));

    double f2 = smooth_mu_f2_x_minus_mu_f1_over_x3(x, mu_s, mu_k, eps_v);
    f2 *= x * x * x;
    f2 += smooth_mu_f1(x, mu_s, mu_k, eps_v);
    f2 /= x;

    CHECK(f2 == Catch::Approx(fd_f2[0]).margin(MARGIN).epsilon(EPSILON));
}