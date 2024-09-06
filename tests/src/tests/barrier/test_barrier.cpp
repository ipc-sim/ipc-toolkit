#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <catch2/generators/catch_generators_random.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>

#include <ipc/barrier/barrier.hpp>

#include <finitediff.hpp>

namespace ipc {

class NormalizedClampedLogBarrier : public ipc::Barrier {
public:
    double operator()(const double d, const double dhat) const override
    {
        if (d <= 0.0) {
            return std::numeric_limits<double>::infinity();
        }
        if (d >= dhat) {
            return 0;
        }

        // b(d) = -(d/d̂-1)²ln(d / d̂)
        const auto d_dhat = d / dhat;
        const auto d_dhat_minus_1 = d_dhat - 1;
        return -d_dhat_minus_1 * d_dhat_minus_1 * std::log(d_dhat);
    }

    double first_derivative(const double d, const double dhat) const override
    {
        if (d <= 0.0 || d >= dhat) {
            return 0.0;
        }
        const double t0 = 1.0 / dhat;
        const double t1 = d * t0;
        const double t2 = 1 - t1;
        return t2 * (2 * t0 * std::log(t1) - t2 / d);
    }

    double second_derivative(const double d, const double dhat) const override
    {
        if (d <= 0.0 || d >= dhat) {
            return 0.0;
        }

        const double t0 = 1.0 / dhat;
        const double t1 = d * t0;
        const double t2 = 1 - t1;
        return 4 * t0 * t2 / d + (t2 * t2) / (d * d)
            - 2 * std::log(t1) / (dhat * dhat);
    }

    double units(const double dhat) const override { return 1; }
};

/// @warning This implementation will not work with dmin > 0
class PhysicalBarrier : public ipc::NormalizedClampedLogBarrier {
public:
    PhysicalBarrier(const bool use_dist_sqr) : use_dist_sqr(use_dist_sqr) { }

    double operator()(const double d, const double dhat) const override
    {
        return (use_dist_sqr ? sqrt(dhat) : dhat)
            * ipc::NormalizedClampedLogBarrier::operator()(d, dhat);
    }

    double first_derivative(const double d, const double dhat) const override
    {
        return (use_dist_sqr ? sqrt(dhat) : dhat)
            * ipc::NormalizedClampedLogBarrier::first_derivative(d, dhat);
    }

    double second_derivative(const double d, const double dhat) const override
    {
        return (use_dist_sqr ? sqrt(dhat) : dhat)
            * ipc::NormalizedClampedLogBarrier::second_derivative(d, dhat);
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