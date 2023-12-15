#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <catch2/generators/catch_generators_random.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>

#include <ipc/barrier/barrier.hpp>

#include <finitediff.hpp>

namespace {
double normalized_barrier(const double d, const double dhat)
{
    if (d <= 0.0) {
        return std::numeric_limits<double>::infinity();
    }
    if (d >= dhat) {
        return 0;
    }

    // b(d) = -(d/d̂-1)²ln(d / d̂)
    const auto t0 = d / dhat;
    return -std::pow(1 - t0, 2) * std::log(t0);
}

double normalized_barrier_first_derivative(const double d, const double dhat)
{
    if (d <= 0.0 || d >= dhat) {
        return 0.0;
    }
    const double t0 = 1.0 / dhat;
    const double t1 = d * t0;
    const double t2 = 1 - t1;
    return t2 * (2 * t0 * std::log(t1) - t2 / d);
}

double normalized_barrier_second_derivative(const double d, const double dhat)
{
    if (d <= 0.0 || d >= dhat) {
        return 0.0;
    }

    const double t0 = 1.0 / dhat;
    const double t1 = d * t0;
    const double t2 = 1 - t1;
    return 4 * t0 * t2 / d + std::pow(t2, 2) / std::pow(d, 2)
        - 2 * std::log(t1) / std::pow(dhat, 2);
}
} // namespace

TEST_CASE("Barrier derivatives", "[barrier]")
{
    bool use_dist_sqr = GENERATE(false, true);
    double dhat = GENERATE_COPY(range(use_dist_sqr ? -2 : -5, 0));
    dhat = pow(10, dhat);

    double d =
        GENERATE_COPY(take(10, random(dhat / 2, 0.9 * dhat))); // ∈ [0, d̂]
    Eigen::Matrix<double, 1, 1> d_vec;
    d_vec << d;

    // Check derivatives

    std::function<double(double, double)> barrier;
    std::function<double(double, double)> barrier_first_derivative;
    std::function<double(double, double)> barrier_second_derivative;
    SECTION("Original IPC barrier")
    {
        barrier = ipc::barrier;
        barrier_first_derivative = ipc::barrier_first_derivative;
        barrier_second_derivative = ipc::barrier_second_derivative;
    }
    SECTION("Normalized barrier")
    {
        barrier = normalized_barrier;
        barrier_first_derivative = normalized_barrier_first_derivative;
        barrier_second_derivative = normalized_barrier_second_derivative;
    }
    SECTION("Barrier with physical units")
    {
        barrier = [dhat](double _d, double p_dhat) {
            return dhat * normalized_barrier(_d, p_dhat);
        };
        barrier_first_derivative = [dhat](double _d, double p_dhat) {
            return dhat * normalized_barrier_first_derivative(_d, p_dhat);
        };
        barrier_second_derivative = [dhat](double _d, double p_dhat) {
            return dhat * normalized_barrier_second_derivative(_d, p_dhat);
        };
    }

    if (use_dist_sqr) {
        d_vec *= d;
        d *= d;
        dhat *= dhat;
    }

    Eigen::VectorXd fgrad(1);
    fd::finite_gradient(
        d_vec, [&](const Eigen::VectorXd& _d) { return barrier(_d[0], dhat); },
        fgrad);

    Eigen::VectorXd grad(1);
    grad << barrier_first_derivative(d, dhat);

    CAPTURE(dhat, d, fgrad(0), grad(0), use_dist_sqr);
    CHECK(fd::compare_gradient(fgrad, grad));

    // Check second_derivative

    fd::finite_gradient(
        d_vec,
        [&](const Eigen::VectorXd& _d) {
            return barrier_first_derivative(_d[0], dhat);
        },
        fgrad);

    grad << barrier_second_derivative(d, dhat);

    CAPTURE(dhat, d, fgrad(0), grad(0), use_dist_sqr);
    CHECK(fd::compare_gradient(fgrad, grad));
}

TEST_CASE("Physical barrier", "[barrier]")
{
    const bool use_dist_sqr = GENERATE(false, true);
    const double dhat =
        pow(10, GENERATE_COPY(range(use_dist_sqr ? -5 : -5, 0)));

    const double d = GENERATE_COPY(take(10, random(dhat / 2, 0.9 * dhat)));

    const double p_d = (use_dist_sqr ? d : 1) * d;
    const double p_dhat = (use_dist_sqr ? dhat : 1) * dhat;
    const double divisor = use_dist_sqr ? std::pow(dhat, 3) : dhat;

    CAPTURE(use_dist_sqr, dhat, d, divisor, p_d, p_dhat);

    const double b_original = ipc::barrier(p_d, p_dhat) / divisor;
    const double b_new = dhat * normalized_barrier(p_d, p_dhat);

    CHECK(b_original == Catch::Approx(b_new));

    const double b_original_first_derivative =
        ipc::barrier_first_derivative(p_d, p_dhat) / divisor;
    const double b_new_first_derivative =
        dhat * normalized_barrier_first_derivative(p_d, p_dhat);

    CHECK(b_original_first_derivative == Catch::Approx(b_new_first_derivative));

    const double b_original_second_derivative =
        ipc::barrier_second_derivative(p_d, p_dhat) / divisor;
    const double b_new_second_derivative =
        dhat * normalized_barrier_second_derivative(p_d, p_dhat);

    CHECK(
        b_original_second_derivative == Catch::Approx(b_new_second_derivative));
}