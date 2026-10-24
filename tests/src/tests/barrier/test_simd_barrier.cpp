#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <tests/simd_utils.hpp>

#ifdef IPC_TOOLKIT_WITH_SIMD

#include <ipc/barrier/barrier.hpp>

#include <array>

using namespace ipc;
using namespace ipc::tests;

namespace {

/// @brief d <= 0 (penetration), d in (0, dhat/2), d in (dhat/2, dhat), and
/// d >= dhat (inactive) -- every branch of every barrier.
std::array<double, 4> regions(const double dhat)
{
    return { { -0.1 * dhat, 0.1 * dhat, 0.7 * dhat, 1.5 * dhat } };
}

/// @brief The barriers take two arguments, so bind dhat and sweep d.
///
/// The tolerance is looser than a value's: a barrier is a log of a ratio, so
/// the two paths' multiply-add contraction shows up amplified.
template <typename BarrierValue, typename BarrierBatch>
void check_matches_scalar(
    const std::string& name,
    const double dhat,
    BarrierValue&& value,
    BarrierBatch&& batch_value)
{
    check_swept_lanes(
        name, regions(dhat), [&](double d) { return value(d, dhat); },
        [&](Batch d) { return batch_value(d, Batch(dhat)); }, 1e-12);
}

/// @brief Compare f(barrier, d, dhat) against f(batch_barrier, d, dhat), where
/// one generic callable serves both the scalar and the batch barrier.
template <typename ScalarBarrier, typename BatchBarrier, typename F>
void check_matches_scalar(
    const std::string& name,
    const double dhat,
    const ScalarBarrier& barrier,
    const BatchBarrier& batch_barrier,
    F&& f)
{
    check_matches_scalar(
        name, dhat, [&](double d, double dh) { return f(barrier, d, dh); },
        [&](Batch d, Batch dh) { return f(batch_barrier, d, dh); });
}

} // namespace

TEST_CASE(
    "SIMD batch barriers match the scalar ones lane-wise", "[barrier][simd]")
{
    const double dhat = 1e-2;

    check_matches_scalar(
        "barrier", dhat, [](double d, double dh) { return barrier(d, dh); },
        [](Batch d, Batch dh) { return barrier(d, dh); });
    check_matches_scalar(
        "barrier_first_derivative", dhat,
        [](double d, double dh) { return barrier_first_derivative(d, dh); },
        [](Batch d, Batch dh) { return barrier_first_derivative(d, dh); });
    check_matches_scalar(
        "barrier_second_derivative", dhat,
        [](double d, double dh) { return barrier_second_derivative(d, dh); },
        [](Batch d, Batch dh) { return barrier_second_derivative(d, dh); });

    const ClampedLogSqBarrier<double> log_sq;
    const ClampedLogSqBarrier<Batch> log_sq_batch;
    check_matches_scalar(
        "ClampedLogSqBarrier", dhat,
        [&](double d, double dh) { return log_sq(d, dh); },
        [&](Batch d, Batch dh) { return log_sq_batch(d, dh); });
    check_matches_scalar(
        "ClampedLogSqBarrier::first_derivative", dhat,
        [&](double d, double dh) { return log_sq.first_derivative(d, dh); },
        [&](Batch d, Batch dh) {
            return log_sq_batch.first_derivative(d, dh);
        });
    check_matches_scalar(
        "ClampedLogSqBarrier::second_derivative", dhat,
        [&](double d, double dh) { return log_sq.second_derivative(d, dh); },
        [&](Batch d, Batch dh) {
            return log_sq_batch.second_derivative(d, dh);
        });

    const CubicBarrier<double> cubic;
    const CubicBarrier<Batch> cubic_batch;
    check_matches_scalar(
        "CubicBarrier", dhat, [&](double d, double dh) { return cubic(d, dh); },
        [&](Batch d, Batch dh) { return cubic_batch(d, dh); });
    check_matches_scalar(
        "CubicBarrier::first_derivative", dhat,
        [&](double d, double dh) { return cubic.first_derivative(d, dh); },
        [&](Batch d, Batch dh) { return cubic_batch.first_derivative(d, dh); });
    check_matches_scalar(
        "CubicBarrier::second_derivative", dhat,
        [&](double d, double dh) { return cubic.second_derivative(d, dh); },
        [&](Batch d, Batch dh) {
            return cubic_batch.second_derivative(d, dh);
        });

    const TwoStageBarrier<double> two_stage;
    const TwoStageBarrier<Batch> two_stage_batch;
    check_matches_scalar(
        "TwoStageBarrier", dhat,
        [&](double d, double dh) { return two_stage(d, dh); },
        [&](Batch d, Batch dh) { return two_stage_batch(d, dh); });
    check_matches_scalar(
        "TwoStageBarrier::first_derivative", dhat,
        [&](double d, double dh) { return two_stage.first_derivative(d, dh); },
        [&](Batch d, Batch dh) {
            return two_stage_batch.first_derivative(d, dh);
        });
    check_matches_scalar(
        "TwoStageBarrier::second_derivative", dhat,
        [&](double d, double dh) { return two_stage.second_derivative(d, dh); },
        [&](Batch d, Batch dh) {
            return two_stage_batch.second_derivative(d, dh);
        });
}

TEST_CASE(
    "SIMD batch inverse-power and near-far barriers match the scalar ones "
    "lane-wise",
    "[barrier][simd]")
{
    const double dhat = 1e-2;
    const double power = GENERATE(1.0, 2.0, 3.0);
    // The far weight ramps up over [alpha dhat / 2, alpha dhat], so 0.1 dhat
    // is on the ramp for alpha = 0.15 and 0.7 dhat is on it for alpha = 1.
    // Between the two, each interior region of regions() is checked both on
    // the ramp and off it.
    const double alpha = GENERATE(0.15, 1.0);
    CAPTURE(power, alpha);

    const auto value = [](const auto& b, auto d, auto dh) { return b(d, dh); };
    const auto first_derivative = [](const auto& b, auto d, auto dh) {
        return b.first_derivative(d, dh);
    };
    const auto second_derivative = [](const auto& b, auto d, auto dh) {
        return b.second_derivative(d, dh);
    };

    const InversePowerBarrier<double> inverse_power(power);
    const InversePowerBarrier<Batch> inverse_power_batch(power);
    check_matches_scalar(
        "InversePowerBarrier", dhat, inverse_power, inverse_power_batch, value);
    check_matches_scalar(
        "InversePowerBarrier::first_derivative", dhat, inverse_power,
        inverse_power_batch, first_derivative);
    check_matches_scalar(
        "InversePowerBarrier::second_derivative", dhat, inverse_power,
        inverse_power_batch, second_derivative);

    // The inverse-power base is +inf at penetration, which is exactly where a
    // batch must not let the blended-away ramp lanes leak a NaN.
    const NearFarBarrier<double> near_far(&inverse_power, alpha);
    const NearFarBarrier<Batch> near_far_batch(&inverse_power_batch, alpha);
    check_matches_scalar(
        "NearFarBarrier::near_value", dhat, near_far, near_far_batch,
        [](const auto& b, auto d, auto dh) { return b.near_value(d, dh); });
    check_matches_scalar(
        "NearFarBarrier::far_value", dhat, near_far, near_far_batch,
        [](const auto& b, auto d, auto dh) { return b.far_value(d, dh); });
    check_matches_scalar(
        "NearFarBarrier::first_derivative_near", dhat, near_far, near_far_batch,
        [](const auto& b, auto d, auto dh) {
            return b.first_derivative_near(d, dh);
        });
    check_matches_scalar(
        "NearFarBarrier::first_derivative_far", dhat, near_far, near_far_batch,
        [](const auto& b, auto d, auto dh) {
            return b.first_derivative_far(d, dh);
        });
    check_matches_scalar(
        "NearFarBarrier::second_derivative_near", dhat, near_far,
        near_far_batch, [](const auto& b, auto d, auto dh) {
            return b.second_derivative_near(d, dh);
        });
    check_matches_scalar(
        "NearFarBarrier::second_derivative_far", dhat, near_far, near_far_batch,
        [](const auto& b, auto d, auto dh) {
            return b.second_derivative_far(d, dh);
        });
}

#endif
