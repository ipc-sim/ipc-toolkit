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

TEST_CASE(
    "anisotropic_mu_eff_f helper",
    "[friction][smooth-mu][anisotropic]")
{
    static constexpr double EPSILON = 1e-6;
    static constexpr double MARGIN = 1e-8;

    // Test with various inputs
    const double mu_s0 = GENERATE(0.1, 0.5, 1.0);
    const double mu_s1 = GENERATE(0.1, 0.5, 1.0);
    const double mu_k0 = GENERATE(0.1, 0.5, 1.0);
    const double mu_k1 = GENERATE(0.1, 0.5, 1.0);

    Eigen::Vector2d mu_s_aniso(mu_s0, mu_s1);
    Eigen::Vector2d mu_k_aniso(mu_k0, mu_k1);

    // Test unit directions
    Eigen::Vector2d tau_dir_x(1.0, 0.0);
    Eigen::Vector2d tau_dir_y(0.0, 1.0);
    Eigen::Vector2d tau_dir_diag(1.0 / std::sqrt(2.0), 1.0 / std::sqrt(2.0));

    // Test x direction
    const auto [mu_s_eff_x, mu_k_eff_x] =
        anisotropic_mu_eff_f(
            tau_dir_x, mu_s_aniso, mu_k_aniso);
    CHECK(
        mu_s_eff_x
        == Catch::Approx(mu_s_aniso[0]).margin(MARGIN).epsilon(EPSILON));
    CHECK(
        mu_k_eff_x
        == Catch::Approx(mu_k_aniso[0]).margin(MARGIN).epsilon(EPSILON));

    // Test y direction
    const auto [mu_s_eff_y, mu_k_eff_y] =
        anisotropic_mu_eff_f(
            tau_dir_y, mu_s_aniso, mu_k_aniso);
    CHECK(
        mu_s_eff_y
        == Catch::Approx(mu_s_aniso[1]).margin(MARGIN).epsilon(EPSILON));
    CHECK(
        mu_k_eff_y
        == Catch::Approx(mu_k_aniso[1]).margin(MARGIN).epsilon(EPSILON));

    // Test diagonal direction
    const auto [mu_s_eff_diag, mu_k_eff_diag] =
        anisotropic_mu_eff_f(
            tau_dir_diag, mu_s_aniso, mu_k_aniso);
    const double expected_mu_s_diag = std::sqrt(
        mu_s_aniso[0] * mu_s_aniso[0] * 0.5
        + mu_s_aniso[1] * mu_s_aniso[1] * 0.5);
    const double expected_mu_k_diag = std::sqrt(
        mu_k_aniso[0] * mu_k_aniso[0] * 0.5
        + mu_k_aniso[1] * mu_k_aniso[1] * 0.5);
    CHECK(
        mu_s_eff_diag
        == Catch::Approx(expected_mu_s_diag).margin(MARGIN).epsilon(EPSILON));
    CHECK(
        mu_k_eff_diag
        == Catch::Approx(expected_mu_k_diag).margin(MARGIN).epsilon(EPSILON));
}

TEST_CASE(
    "anisotropic_mu_eff_f_dtau helper",
    "[friction][smooth-mu][anisotropic-derivative]")
{
    static constexpr double EPSILON = 1e-4;
    static constexpr double MARGIN = 1e-6;
    static constexpr double H = 1e-8;

    const double mu0 = 0.5;
    const double mu1 = 0.8;
    Eigen::Vector2d mu_aniso(mu0, mu1);

    // Test with non-zero tau
    Eigen::Vector2d tau(0.5, 0.3);
    constexpr double tiny = 1e-10;
    Eigen::Vector2d tau_dir;
    if (tau.norm() < tiny) {
        tau_dir = Eigen::Vector2d(1.0, 0.0);
    } else {
        tau_dir = tau / tau.norm();
    }
    const auto [mu_s_eff, mu_k_eff] =
        anisotropic_mu_eff_f(
            tau_dir, mu_aniso, mu_aniso);
    const double mu_eff = mu_s_eff;

    const Eigen::Vector2d dmu_eff_dtau =
        anisotropic_mu_eff_f_dtau(tau, mu_aniso, mu_eff);

    // Finite difference check
    Eigen::Matrix<double, 2, 1> Tau;
    Tau << tau[0], tau[1];

    Eigen::VectorXd fd_dmu_eff_dtau(2);
    fd::finite_gradient(
        Tau,
        [&](const Eigen::VectorXd& _Tau) {
            Eigen::Vector2d _tau(_Tau[0], _Tau[1]);
            Eigen::Vector2d _tau_dir;
            if (_tau.norm() < tiny) {
                _tau_dir = Eigen::Vector2d(1.0, 0.0);
            } else {
                _tau_dir = _tau / _tau.norm();
            }
            const auto [_mu_s_eff, _mu_k_eff] =
                anisotropic_mu_eff_f(
                    _tau_dir, mu_aniso, mu_aniso);
            return _mu_s_eff;
        },
        fd_dmu_eff_dtau, fd::AccuracyOrder::SECOND, H);

    CHECK(
        dmu_eff_dtau[0]
        == Catch::Approx(fd_dmu_eff_dtau[0]).margin(MARGIN).epsilon(EPSILON));
    CHECK(
        dmu_eff_dtau[1]
        == Catch::Approx(fd_dmu_eff_dtau[1]).margin(MARGIN).epsilon(EPSILON));

    // Test edge case: ||tau|| ≈ 0
    Eigen::Vector2d tau_zero(1e-12, 1e-12);
    Eigen::Vector2d tau_dir_zero = tau_zero / tau_zero.norm();
    const auto [mu_s_eff_zero, mu_k_eff_zero] =
        anisotropic_mu_eff_f(
            tau_dir_zero, mu_aniso, mu_aniso);
    const Eigen::Vector2d dmu_eff_dtau_zero =
        anisotropic_mu_eff_f_dtau(tau_zero, mu_aniso, mu_s_eff_zero);
    CHECK(dmu_eff_dtau_zero.norm() < 1e-6); // Should be approximately zero
}