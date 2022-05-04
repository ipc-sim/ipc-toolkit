#include <Eigen/Core>
#include <catch2/catch.hpp>

#include <finitediff.hpp>

#include <ipc/barrier/barrier.hpp>

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
        barrier = ipc::physical_barrier<double>;
        barrier_gradient = ipc::physical_barrier_gradient;
        barrier_hessian = ipc::physical_barrier_hessian;
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
