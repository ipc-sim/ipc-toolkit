#include <catch2/catch_test_macros.hpp>

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
    return { -0.1 * dhat, 0.1 * dhat, 0.7 * dhat, 1.5 * dhat };
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

#endif
