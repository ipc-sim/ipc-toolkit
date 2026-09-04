#include <Eigen/Core>
#include <catch2/catch_all.hpp>

#include <ipc/dynamics/time_integration/bdf.hpp>

using namespace ipc;
using namespace ipc::dynamics;

TEST_CASE("BDF coefficients", "[dynamics][bdf]")
{
    // The α coefficients of each order sum to 1 (consistency).
    for (int order = 1; order <= 6; ++order) {
        const std::vector<double>& alpha = BDF::alphas(order - 1);
        REQUIRE(alpha.size() == order);
        double sum = 0;
        for (const double a : alpha) {
            sum += a;
        }
        CHECK(sum == Catch::Approx(1.0));
        CHECK(BDF::betas(order - 1) > 0.0);
    }
}

TEST_CASE("BDF order ramps up", "[dynamics][bdf]")
{
    const int target = GENERATE(1, 2, 3, 4, 5, 6);

    BDF bdf(target);
    bdf.set_dt(0.01);
    const int n = 3; // DOFs (single "body" worth is irrelevant here)
    bdf.init(
        Eigen::VectorXd::Zero(n), Eigen::VectorXd::Zero(n),
        Eigen::VectorXd::Zero(n), /*num_bodies=*/1);

    CHECK(bdf.target_order() == target);

    // The effective order starts at 1 and increases by one per accepted step
    // until it reaches the target (BDF startup).
    CHECK(bdf.order() == 1);
    for (int step = 1; step <= target + 2; ++step) {
        bdf.update(Eigen::VectorXd::Zero(n));
        CHECK(bdf.order() == std::min(target, step + 1));
    }
}

TEST_CASE("BDF is exact for constant velocity", "[dynamics][bdf]")
{
    // BDF-n reproduces polynomials of degree ≤ n exactly. For a
    // constant-velocity trajectory x(t) = x0 + v0 t, once enough history has
    // accumulated compute_velocity() must recover v0 and predicted_positions()
    // must extrapolate exactly — for any order.
    const int target = GENERATE(1, 2, 3, 4);
    const int n = 5;
    const double dt = 0.01;

    const Eigen::VectorXd x0 = Eigen::VectorXd::Random(n);
    const Eigen::VectorXd v0 = Eigen::VectorXd::Random(n);

    BDF bdf(target);
    bdf.set_dt(dt);
    bdf.init(x0, v0, Eigen::VectorXd::Zero(n), /*num_bodies=*/1);

    Eigen::VectorXd x = x0;
    // Warm up the history along the exact linear trajectory.
    for (int step = 1; step <= target; ++step) {
        x = x0 + step * dt * v0;
        bdf.update(x);
    }

    // Now the effective order equals the target; the scheme is exact.
    REQUIRE(bdf.order() == target);

    const Eigen::VectorXd x_next = x + dt * v0;
    CHECK((bdf.compute_velocity(x_next) - v0).norm() < 1e-10);
    CHECK((bdf.predicted_positions() - x_next).norm() < 1e-10);

    // acceleration_scaling() = (β Δt)²
    const double beta_dt = BDF::betas(target - 1) * dt;
    CHECK(bdf.acceleration_scaling() == Catch::Approx(beta_dt * beta_dt));
    CHECK(bdf.velocity_scaling() == Catch::Approx(beta_dt));
}

TEST_CASE("BDF invalid order throws", "[dynamics][bdf]")
{
    CHECK_THROWS(BDF(0));
    CHECK_THROWS(BDF(7));
    CHECK_NOTHROW(BDF(1));
    CHECK_NOTHROW(BDF(6));
}
