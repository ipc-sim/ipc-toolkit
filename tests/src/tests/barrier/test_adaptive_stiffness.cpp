#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <ipc/barrier/adaptive_stiffness.hpp>

using namespace ipc;

TEST_CASE("Initial barrier stiffness", "[adaptive_stiffness]")
{
    const double bbox_diagonal = 1.0;
    const ClampedLogBarrier barrier;
    const double dhat = 1e-3;
    const double average_mass = 1.0;
    const Eigen::VectorXd grad_energy = Eigen::VectorXd::Constant(1, 100);
    const Eigen::VectorXd grad_barrier = Eigen::VectorXd::Constant(1, -100);

    double max_barrier_stiffness;
    const double kappa = initial_barrier_stiffness(
        bbox_diagonal, barrier, dhat, average_mass, grad_energy, grad_barrier,
        max_barrier_stiffness);

    const double expected_kappa = 1e11 * average_mass
        / (4 * 1e-16 * barrier_second_derivative(1e-16, dhat * dhat));

    CHECK(kappa == Catch::Approx(expected_kappa));
    CHECK(max_barrier_stiffness == Catch::Approx(100 * expected_kappa));
}

TEST_CASE("Update barrier stiffness", "[adaptive_stiffness]")
{
    double prev_min_distance = 1e-10 * 1e-10;
    double min_distance = 1e-11 * 1e-11;
    double max_barrier_stiffness = 10.0;
    double barrier_stiffness = 1.0;
    double bbox_diagonal = 1.0;

    double kappa = update_barrier_stiffness(
        prev_min_distance, min_distance, max_barrier_stiffness,
        barrier_stiffness, bbox_diagonal);

    REQUIRE(kappa == Catch::Approx(2.0));
}