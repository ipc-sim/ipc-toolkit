#include <Eigen/Core>
#include <catch2/catch.hpp>

#include <finitediff.hpp>

#include <ipc/barrier/barrier.hpp>

namespace {
double physical_barrier(const double d, const double dhat)
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
    const double d_over_dhat = d / dhat;
    const double d_over_dhat_minus_1 = d_over_dhat - 1;
    return -dhat * d_over_dhat_minus_1 * d_over_dhat_minus_1 * log(d_over_dhat);
}

double physical_barrier_gradient(const double d, const double dhat)
{
    if (d <= 0.0 || d >= dhat) {
        return 0.0;
    }
    return (d - dhat) * (-2 * d * log(d / dhat) - d + dhat) / (d * dhat);
}

double physical_barrier_hessian(const double d, const double dhat)
{
    if (d <= 0.0 || d >= dhat) {
        return 0.0;
    }
    const double d_minus_dhat = d - dhat;
    return (-2 * log(d / dhat) + (d_minus_dhat * d_minus_dhat) / (d * d) - 4)
        / dhat
        + 4 / d;
}
} // namespace

TEST_CASE("Test barrier derivatives", "[barrier]")
{
    double dhat = GENERATE(range(-5, 2));
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
        barrier = ipc::barrier<double>;
        barrier_gradient = ipc::barrier_gradient;
        barrier_hessian = ipc::barrier_hessian;
    }
    SECTION("Barrier with physical units")
    {
        barrier = physical_barrier;
        barrier_gradient = physical_barrier_gradient;
        barrier_hessian = physical_barrier_hessian;
    }

    Eigen::VectorXd fgrad(1);
    fd::finite_gradient(
        d_vec, [&](const Eigen::VectorXd& d) { return barrier(d[0], dhat); },
        fgrad);

    Eigen::VectorXd grad(1);
    grad << barrier_gradient(d, dhat);

    CAPTURE(dhat, d, fgrad(0), grad(0));
    CHECK(fd::compare_gradient(fgrad, grad));

    // Check hessian

    fd::finite_gradient(
        d_vec,
        [&](const Eigen::VectorXd& d) { return barrier_gradient(d[0], dhat); },
        fgrad);

    grad << barrier_hessian(d, dhat);

    CAPTURE(dhat, d, fgrad(0), grad(0));
    CHECK(fd::compare_gradient(fgrad, grad));
}

TEST_CASE("Test physical barrier", "[barrier]")
{
    double dhat = GENERATE(range(-5, 2));
    dhat = pow(10, dhat);

    double d =
        GENERATE_COPY(take(10, random(dhat / 2, 0.9 * dhat))); // ∈ [0, d̂]

    double b_original = ipc::barrier(d, dhat) / dhat;
    double b_new = physical_barrier(d, dhat);

    CHECK(b_original == Approx(b_new));

    double b_original_gradient = ipc::barrier_gradient(d, dhat) / dhat;
    double b_new_gradient = physical_barrier_gradient(d, dhat);

    CHECK(b_original_gradient == Approx(b_new_gradient));

    double b_original_hessian = ipc::barrier_hessian(d, dhat) / dhat;
    double b_new_hessian = physical_barrier_hessian(d, dhat);

    CHECK(b_original_hessian == Approx(b_new_hessian));
}