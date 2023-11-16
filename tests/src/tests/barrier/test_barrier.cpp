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

double normalized_barrier_gradient(const double d, const double dhat)
{
    if (d <= 0.0 || d >= dhat) {
        return 0.0;
    }
    const double t0 = 1.0 / dhat;
    const double t1 = d * t0;
    const double t2 = 1 - t1;
    return t2 * (2 * t0 * std::log(t1) - t2 / d);
}

double normalized_barrier_hessian(const double d, const double dhat)
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

    // Check gradient

    std::function<double(double, double)> barrier, barrier_gradient,
        barrier_hessian;
    SECTION("Original IPC barrier")
    {
        barrier = ipc::barrier;
        barrier_gradient = ipc::barrier_gradient;
        barrier_hessian = ipc::barrier_hessian;
    }
    SECTION("Normalized barrier")
    {
        barrier = normalized_barrier;
        barrier_gradient = normalized_barrier_gradient;
        barrier_hessian = normalized_barrier_hessian;
    }
    SECTION("Barrier with physical units")
    {
        barrier = [dhat](double _d, double p_dhat) {
            return dhat * normalized_barrier(_d, p_dhat);
        };
        barrier_gradient = [dhat](double _d, double p_dhat) {
            return dhat * normalized_barrier_gradient(_d, p_dhat);
        };
        barrier_hessian = [dhat](double _d, double p_dhat) {
            return dhat * normalized_barrier_hessian(_d, p_dhat);
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
    grad << barrier_gradient(d, dhat);

    CAPTURE(dhat, d, fgrad(0), grad(0), use_dist_sqr);
    CHECK(fd::compare_gradient(fgrad, grad));

    // Check hessian

    fd::finite_gradient(
        d_vec,
        [&](const Eigen::VectorXd& _d) {
            return barrier_gradient(_d[0], dhat);
        },
        fgrad);

    grad << barrier_hessian(d, dhat);

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

    const double b_original_gradient =
        ipc::barrier_gradient(p_d, p_dhat) / divisor;
    const double b_new_gradient =
        dhat * normalized_barrier_gradient(p_d, p_dhat);

    CHECK(b_original_gradient == Catch::Approx(b_new_gradient));

    const double b_original_hessian =
        ipc::barrier_hessian(p_d, p_dhat) / divisor;
    const double b_new_hessian = dhat * normalized_barrier_hessian(p_d, p_dhat);

    CHECK(b_original_hessian == Catch::Approx(b_new_hessian));
}