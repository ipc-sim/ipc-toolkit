#include <catch2/catch_test_macros.hpp>

#include <tests/simd_utils.hpp>

#ifdef IPC_TOOLKIT_WITH_SIMD

#include <cmath>
#include <limits>

using namespace ipc::tests;

// Every SIMD test compares through `approx`, so its contract is worth pinning
// down here rather than rediscovering it as a mysteriously passing test.

TEST_CASE("SIMD test approx zeroes Approx's default epsilon", "[utils][simd]")
{
    // The reason `approx` sets epsilon(0): Approx ORs its margin test with an
    // epsilon test whose default is float-grade (FLT_EPSILON * 100, about
    // 1.19e-5 relative). Left in place, it would accept a value that misses by
    // far more than VALUE_TOL, which is exactly the bug this guards.
    const double x = 1.0;
    const double off_by_1e7 = x + 1e-7;

    CHECK(off_by_1e7 != approx(x, VALUE_TOL));
    CHECK(off_by_1e7 == Catch::Approx(x)); // the default that must not leak in
}

TEST_CASE(
    "SIMD test approx bounds relative to the given scale", "[utils][simd]")
{
    // The scale is the caller's choice -- for a derivative it is the magnitude
    // of the whole result, so an entry that is a cancellation of much larger
    // terms is not held to its own tiny magnitude.
    CHECK(1e-8 == approx(0.0, DERIVATIVE_TOL, /*scale=*/1e3));
    CHECK(1e-8 != approx(0.0, DERIVATIVE_TOL, /*scale=*/1.0));
}

TEST_CASE("SIMD test approx matches infinities by sign", "[utils][simd]")
{
    // Approx compares without subtracting, which is what makes this work. The
    // barrier tests depend on it: a barrier at penetration must resolve to
    // +inf or 0, never NaN.
    constexpr double INF = std::numeric_limits<double>::infinity();
    constexpr double NAN_ = std::numeric_limits<double>::quiet_NaN();

    CHECK(INF == approx(INF, VALUE_TOL));
    CHECK(-INF == approx(-INF, VALUE_TOL));
    CHECK(-INF != approx(INF, VALUE_TOL));
    CHECK(1.0 != approx(INF, VALUE_TOL));
    CHECK(INF != approx(1.0, VALUE_TOL));

    // A NaN never passes, on either side.
    CHECK(NAN_ != approx(INF, VALUE_TOL));
    CHECK(NAN_ != approx(1.0, VALUE_TOL));
    CHECK(1.0 != approx(NAN_, VALUE_TOL));
}

TEST_CASE("SIMD test Lanes round-trip through a batch", "[utils][simd]")
{
    Lanes in;
    for (int l = 0; l < L; ++l) {
        in[l] = 1.0 + 0.25 * l;
    }

    const Lanes out = Lanes::from(in.load());
    for (int l = 0; l < L; ++l) {
        CAPTURE(l);
        CHECK(out[l] == in[l]);
    }
}

TEST_CASE("SIMD test lane_cases wraps round-robin", "[utils][simd]")
{
    // Two cases and (usually) more lanes, so the wrap is exercised whatever
    // the batch width; the offset rotates which case each lane gets.
    const std::array<double, 2> cases = { { 10.0, 20.0 } };

    for (int offset = 0; offset < 2; ++offset) {
        CAPTURE(offset);
        const Lanes lanes = lane_cases(cases, offset);
        for (int l = 0; l < L; ++l) {
            CAPTURE(l);
            CHECK(lanes[l] == cases[size_t(l + offset) % cases.size()]);
        }
    }
}

#endif
