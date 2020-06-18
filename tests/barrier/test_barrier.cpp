#include <Eigen/Core>
#include <catch2/catch.hpp>

#include <finitediff.hpp>

#include <barrier/barrier.hpp>

TEST_CASE("Test barrier derivatives", "[barrier]")
{
    double dhat = GENERATE(range(-5, 2));
    dhat = pow(10, dhat);

    double d =
        GENERATE_COPY(take(10, random(dhat / 2, 0.9 * dhat))); // ∈ [0, d̂]
    Eigen::Matrix<double, 1, 1> d_vec;
    d_vec << d;

    // Check gradient

    Eigen::VectorXd fgrad(1);
    fd::finite_gradient(
        d_vec,
        [&](const Eigen::VectorXd& d) { return ipc::barrier(d[0], dhat); },
        fgrad);

    Eigen::VectorXd grad(1);
    grad << ipc::barrier_gradient(d, dhat);

    CAPTURE(dhat, d, fgrad(0), grad(0));
    CHECK(fd::compare_gradient(fgrad, grad));

    // Check hessian

    fd::finite_gradient(
        d_vec,
        [&](const Eigen::VectorXd& d) {
            return ipc::barrier_gradient(d[0], dhat);
        },
        fgrad);

    grad << ipc::barrier_hessian(d, dhat);

    CAPTURE(dhat, d, fgrad(0), grad(0));
    CHECK(fd::compare_gradient(fgrad, grad));
}
