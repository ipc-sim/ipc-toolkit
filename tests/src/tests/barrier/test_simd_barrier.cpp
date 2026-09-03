#include <catch2/catch_test_macros.hpp>

#include <ipc/utils/simd.hpp>

#ifdef IPC_TOOLKIT_WITH_SIMD

#include <ipc/barrier/barrier.hpp>

#include <cmath>
#include <limits>
#include <string>
#include <vector>

using namespace ipc;

namespace {

using Batch = SimdBatch<double>;
constexpr int L = int(Batch::size);

/// @brief d <= 0 (penetration), d in (0, dhat/2), d in (dhat/2, dhat), and
/// d >= dhat (inactive) — every branch of every barrier.
std::vector<double> regions(const double dhat)
{
    return { -0.1 * dhat, 0.1 * dhat, 0.7 * dhat, 1.5 * dhat };
}

/// @brief Fill the lanes from `regions` starting at `offset`. A batch may hold
/// fewer lanes than there are regions, so the caller sweeps the offset to reach
/// every region.
std::vector<double> lane_ds(const double dhat, const int offset)
{
    const std::vector<double> region_ds = regions(dhat);
    std::vector<double> ds(L);
    for (int l = 0; l < L; ++l) {
        ds[l] = region_ds[(l + offset) % region_ds.size()];
    }
    return ds;
}

bool close(const double got, const double want)
{
    // NaN is deliberately not tolerated: every barrier must resolve d <= 0 to
    // +inf (value) or 0 (derivative), never NaN.
    if (std::isinf(want)) {
        return std::isinf(got) && (got > 0) == (want > 0);
    }
    return std::abs(got - want) <= 1e-12 * std::max(1.0, std::abs(want));
}

template <typename BarrierValue, typename BarrierBatch>
void check_matches_scalar(
    const std::string& name,
    const double dhat,
    BarrierValue&& value,
    BarrierBatch&& batch_value)
{
    for (int offset = 0; offset < int(regions(dhat).size()); ++offset) {
        const std::vector<double> ds = lane_ds(dhat, offset);

        const Batch d_batch = Batch::load_unaligned(ds.data());
        const Batch dhat_batch(dhat);
        const Batch got_batch = batch_value(d_batch, dhat_batch);

        std::vector<double> got(L);
        got_batch.store_unaligned(got.data());

        for (int l = 0; l < L; ++l) {
            const double want = value(ds[l], dhat);
            CAPTURE(name, offset, l, ds[l], dhat, got[l], want);
            CHECK(close(got[l], want));
        }
    }
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
