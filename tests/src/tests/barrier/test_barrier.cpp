#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <catch2/generators/catch_generators_random.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>

#include <ipc/barrier/barrier.hpp>

#include <finitediff.hpp>

#include <memory>

namespace ipc {

/// @warning This implementation will not work with dmin > 0
class PhysicalBarrier : public NormalizedClampedLogBarrier {
public:
    PhysicalBarrier(const bool _use_dist_sqr) : use_dist_sqr(_use_dist_sqr) { }

    double operator()(const double d, const double dhat) const override
    {
        return (use_dist_sqr ? sqrt(dhat) : dhat)
            * NormalizedClampedLogBarrier::operator()(d, dhat);
    }

    double first_derivative(const double d, const double dhat) const override
    {
        return (use_dist_sqr ? sqrt(dhat) : dhat)
            * NormalizedClampedLogBarrier::first_derivative(d, dhat);
    }

    double second_derivative(const double d, const double dhat) const override
    {
        return (use_dist_sqr ? sqrt(dhat) : dhat)
            * NormalizedClampedLogBarrier::second_derivative(d, dhat);
    }

private:
    bool use_dist_sqr;
};

} // namespace ipc

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

    std::unique_ptr<ipc::Barrier> barrier;
    SECTION("Original")
    {
        barrier = std::make_unique<ipc::ClampedLogBarrier>();
    }
    SECTION("Normalized")
    {
        barrier = std::make_unique<ipc::NormalizedClampedLogBarrier>();
    }
    SECTION("Physical")
    {
        barrier = std::make_unique<ipc::PhysicalBarrier>(use_dist_sqr);
    }
    SECTION("ClampedLogSq")
    {
        barrier = std::make_unique<ipc::ClampedLogSqBarrier>();
    }
    SECTION("Cubic") { barrier = std::make_unique<ipc::CubicBarrier>(); }

    if (use_dist_sqr) {
        d_vec *= d;
        d *= d;
        dhat *= dhat;
    }

    Eigen::VectorXd fgrad(1);
    fd::finite_gradient(
        d_vec,
        [&](const Eigen::VectorXd& _d) { return (*barrier)(_d[0], dhat); },
        fgrad);

    Eigen::VectorXd grad(1);
    grad << barrier->first_derivative(d, dhat);

    CAPTURE(dhat, d, fgrad(0), grad(0), use_dist_sqr);
    CHECK(fd::compare_gradient(fgrad, grad));

    // Check second_derivative

    fd::finite_gradient(
        d_vec,
        [&](const Eigen::VectorXd& _d) {
            return barrier->first_derivative(_d[0], dhat);
        },
        fgrad);

    grad << barrier->second_derivative(d, dhat);

    CAPTURE(dhat, d, fgrad(0), grad(0), use_dist_sqr);
    CHECK(fd::compare_gradient(fgrad, grad));
}

TEST_CASE("Physical barrier", "[barrier]")
{
    const bool use_dist_sqr = GENERATE(false, true);

    ipc::ClampedLogBarrier original_barrier;
    ipc::PhysicalBarrier new_barrier(use_dist_sqr);

    const double dhat =
        pow(10, GENERATE_COPY(range(use_dist_sqr ? -5 : -5, 0)));

    const double d = GENERATE_COPY(take(10, random(dhat / 2, 0.9 * dhat)));

    const double p_d = (use_dist_sqr ? d : 1) * d;
    const double p_dhat = (use_dist_sqr ? dhat : 1) * dhat;
    const double divisor = original_barrier.units(p_dhat) / dhat;

    CAPTURE(use_dist_sqr, dhat, d, divisor, p_d, p_dhat);

    const double b_original = original_barrier(p_d, p_dhat) / divisor;
    const double b_new = new_barrier(p_d, p_dhat);

    CHECK(b_original == Catch::Approx(b_new));

    const double b_original_first_derivative =
        original_barrier.first_derivative(p_d, p_dhat) / divisor;
    const double b_new_first_derivative =
        new_barrier.first_derivative(p_d, p_dhat);

    CHECK(b_original_first_derivative == Catch::Approx(b_new_first_derivative));

    const double b_original_second_derivative =
        original_barrier.second_derivative(p_d, p_dhat) / divisor;
    const double b_new_second_derivative =
        new_barrier.second_derivative(p_d, p_dhat);

    CHECK(
        b_original_second_derivative == Catch::Approx(b_new_second_derivative));
}