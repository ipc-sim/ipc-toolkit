#include <Eigen/Core>
#include <catch2/catch_all.hpp>

#include <finitediff.hpp>

#include <ipc/barrier/barrier.hpp>

namespace {
double normalized_barrier(const double d, const double dhat)
{
    // units(d) = m and units(d̂) = m ⟹ units(b(d)) = m
    // units(κ) = Pa ⟹ units(κ b(d)) = Pa m

    if (d <= 0.0) {
        return std::numeric_limits<double>::infinity();
    }
    if (d >= dhat) {
        return 0;
    }

    // b(d) = -d̂(d/d̂-1)²ln(d / d̂)
    const double t0 = d / dhat;
    return -dhat * std::pow(1 - t0, 2) * std::log(t0);
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

TEST_CASE("Test barrier derivatives", "[barrier]")
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
    // SECTION("Normalized barrier")
    // {
    //     barrier = normalized_barrier;
    //     barrier_gradient = normalized_barrier_gradient;
    //     barrier_hessian = normalized_barrier_hessian;
    // }
    SECTION("Barrier with physical units")
    {
        barrier = [dhat](double d, double p_dhat) {
            return dhat * normalized_barrier(d, p_dhat);
        };
        barrier_gradient = [dhat](double d, double p_dhat) {
            return dhat * normalized_barrier_gradient(d, p_dhat);
        };
        barrier_hessian = [dhat](double d, double p_dhat) {
            return dhat * normalized_barrier_hessian(d, p_dhat);
        };
    }

    if (use_dist_sqr) {
        d_vec *= d;
        d *= d;
        dhat *= dhat;
    }

    Eigen::VectorXd fgrad(1);
    fd::finite_gradient(
        d_vec, [&](const Eigen::VectorXd& d) { return barrier(d[0], dhat); },
        fgrad);

    Eigen::VectorXd grad(1);
    grad << barrier_gradient(d, dhat);

    CAPTURE(dhat, d, fgrad(0), grad(0), use_dist_sqr);
    CHECK(fd::compare_gradient(fgrad, grad));

    // Check hessian

    fd::finite_gradient(
        d_vec,
        [&](const Eigen::VectorXd& d) { return barrier_gradient(d[0], dhat); },
        fgrad);

    grad << barrier_hessian(d, dhat);

    CAPTURE(dhat, d, fgrad(0), grad(0), use_dist_sqr);
    CHECK(fd::compare_gradient(fgrad, grad));
}

TEST_CASE("Test physical barrier", "[barrier]")
{
    double dhat = GENERATE(range(-5, 2));
    dhat = pow(10, dhat);

    double d =
        GENERATE_COPY(take(10, random(dhat / 2, 0.9 * dhat))); // ∈ [0, d̂]

    double b_original = ipc::barrier(d, dhat) / dhat;
    double b_new = dhat * normalized_barrier(d, dhat);

    CHECK(b_original == Catch::Approx(b_new));

    double b_original_gradient = ipc::barrier_gradient(d, dhat) / dhat;
    double b_new_gradient = dhat * normalized_barrier_gradient(d, dhat);

    CHECK(b_original_gradient == Catch::Approx(b_new_gradient));

    double b_original_hessian = ipc::barrier_hessian(d, dhat) / dhat;
    double b_new_hessian = dhat * normalized_barrier_hessian(d, dhat);

    CHECK(b_original_hessian == Catch::Approx(b_new_hessian));
}