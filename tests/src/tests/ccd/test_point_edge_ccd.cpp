#include <tests/ccd/collision_generator.hpp>

#include <catch2/catch_test_macros.hpp>

#include <ipc/config.hpp>
#include <ipc/ccd/ccd.hpp>
#include <ipc/ccd/additive_ccd.hpp>

using namespace ipc;

namespace {
/// Compares the time of impact of different implementations
/// against the expected time of impact
void check_toi(
    const VectorMax3d& p_t0,
    const VectorMax3d& e0_t0,
    const VectorMax3d& e1_t0,
    const VectorMax3d& p_t1,
    const VectorMax3d& e0_t1,
    const VectorMax3d& e1_t1,
    const double toi_expected)
{
    double toi;
    bool is_colliding = point_edge_ccd(
        p_t0, e0_t0, e1_t0, p_t1, e0_t1, e1_t1, toi,
        /*min_distance=*/0.0, /*tmax=*/1.0, DEFAULT_CCD_TOLERANCE,
        DEFAULT_CCD_MAX_ITERATIONS, /*conservative_rescaling=*/1.0);
    CHECK(is_colliding);
    CHECK(toi <= toi_expected);

    is_colliding = additive_ccd::point_edge_ccd(
        p_t0, e0_t0, e1_t0, p_t1, e0_t1, e1_t1, toi,
        /*min_distance=*/0.0, /*tmax=*/1.0, /*conservative_rescaling=*/0.99);
    CHECK(is_colliding);
    CHECK(toi <= toi_expected);
}
} // namespace

TEST_CASE("Point-edge CCD", "[ccd][point-edge]")
{
    int dim = GENERATE(2, 3);
    VectorMax3d p_t0(2), p_t1(2), e0_t0(2);
    VectorMax3d e0_t1(2), e1_t0(2), e1_t1(2);
    double toi_expected;
    bool is_collision_expected;
    SECTION("Edge becomes degenerate")
    {
        p_t0 << 0, 1;
        e0_t0 << -1, 0;
        e1_t0 << 1, 0;
        SECTION("Edge degenerates before impact")
        {
            p_t1 << 0, -1;
            e0_t1 << 3, 0;
            e1_t1 << -3, 0;
            // The edge will become degenerate at t=0.25
        }
        SECTION("Edge degenerates after impact")
        {
            p_t1 << 0, -1;
            e0_t1 << 0.5, 0;
            e1_t1 << -0.5, 0;
            // The edge will become degenerate at t=2/3
        }
        // The point will collide with the edge at t=0.5
        is_collision_expected = true;
        toi_expected = 0.5;
    }
    SECTION("Edge moving right; point moving left")
    {
        p_t0 << -1, 0;
        e0_t0 << 1, -1;
        e1_t0 << 1, 1;

        p_t1 << 1, 0;
        e0_t1 << -1, -1;
        e1_t1 << -1, 1;

        is_collision_expected = true;
        toi_expected = 0.5;
    }
    SECTION("Point on edge's line moving towards edge")
    {
        p_t0 << 0, 0;
        e0_t0 << 0, 1;
        e1_t0 << 0, 2;

        p_t1 << 0, 2;
        e0_t1 << 0, 1;
        e1_t1 << 0, 2;

        is_collision_expected = true;
        toi_expected = 0.5;
    }
    SECTION("Point and edge moving parallel")
    {
        p_t0 << 0, 1;
        e0_t0 << 1, 0;
        e1_t0 << 1, 2;

        p_t1 << 0, 2;
        e0_t1 << 1, 1;
        e1_t1 << 1, 3;

        is_collision_expected = false;
    }
    SECTION("Point moving right; edge stretching vertically")
    {
        e1_t0 << 1, -1;
        e1_t1 << 1, -2;
        SECTION("Swap vertices order e_0 = [0, 2]")
        {
            p_t0 << 1, 1;
            e0_t0 << 0, 0;

            p_t1 << 1, 2;
            e0_t1 << 1, 0;

            std::swap(p_t0, e0_t0);
            std::swap(p_t1, e0_t1);
        }
        SECTION("Swap vertices order e_0 = [1, 2]")
        {
            p_t0 << 0, 0;
            e0_t0 << 1, 1;

            p_t1 << 1, 0;
            e0_t1 << 1, 2;
        }
        is_collision_expected = true;
        toi_expected = 1.0;
    }
    SECTION("Point-point")
    {
        p_t0 << 1.11111, 0.5;
        p_t1 << 0.888889, 0.5;
        e0_t0 << 1, 0.5;
        e1_t0 << 1, 0.75;
        e0_t1 = e0_t0;
        e1_t1 = e1_t0;

        is_collision_expected = true;
        toi_expected = 0.5;
    }

    if (dim == 3) {
        p_t0.conservativeResize(3);
        p_t1.conservativeResize(3);
        e0_t0.conservativeResize(3);
        e0_t1.conservativeResize(3);
        e1_t0.conservativeResize(3);
        e1_t1.conservativeResize(3);
        p_t0[2] = p_t1[2] = e0_t0[2] = e1_t0[2] = e0_t1[2] = e1_t1[2] = 0;
    }

    CAPTURE(
        p_t0.transpose(), e0_t0.transpose(), e1_t0.transpose(),
        p_t1.transpose(), e0_t1.transpose(), e1_t1.transpose());

    double toi;
    bool is_colliding;
    // = point_edge_ccd(
    //     p_t0, e0_t0, e1_t0, p_t1, e0_t1, e1_t1, toi,
    //     /*min_distance=*/0.0,
    //     /*tmax=*/1.0, DEFAULT_CCD_TOLERANCE, DEFAULT_CCD_MAX_ITERATIONS,
    //     /*conservative_rescaling=*/1.0);
    // REQUIRE(is_colliding == is_collision_expected);
    // if (is_collision_expected) {
    //     CHECK(toi <= toi_expected);
    // }

    is_colliding = additive_ccd::point_edge_ccd(
        p_t0, e0_t0, e1_t0, p_t1, e0_t1, e1_t1, toi,
        /*min_distance=*/0.0, /*tmax=*/1.0, /*conservative_rescaling=*/0.99);
    REQUIRE(is_colliding == is_collision_expected);
    if (is_collision_expected) {
        CHECK(toi <= toi_expected);
    }
}

#ifdef IPC_TOOLKIT_WITH_INEXACT_CCD
TEST_CASE("Point-edge ToI", "[ccd][point-edge][toi][!mayfail]")
#else
TEST_CASE("Point-edge ToI", "[ccd][point-edge][toi]")
#endif
{
    const int dim = GENERATE(2, 3);
    VectorMax3d p_t0(2), e0_t0(2), e1_t0(2);
    VectorMax3d dp(2), de0(2), de1(2);

    double expected_toi = -1;

    SECTION("perpendicular") // (alpha=0.5)
    {
        p_t0 << 0.0, 1.0;
        e0_t0 << -1.0, 0.0;
        e1_t0 << 1.0, 0.0;

        // touches, intersects, passes-trough
        auto vel = GENERATE(1.0 + 1e-6, 2.0, 4.0);
        // moving direction:
        //  -> ->, -> -, -> <-, - <-, <- <-
        auto j = GENERATE(0, 1, 2, 3, 4);
        // extension, no-deform, compression,
        auto dx = GENERATE(0.5, 0.0, -0.5);

        dp << 0.0, -(3 - j) * vel / 2.0;
        de0 << -dx, (j - 1.0) * vel / 2.0;
        de1 << dx, (j - 1.0) * vel / 2.0;

        expected_toi = 1.0 / vel;
    }

    SECTION("double impact") // (rotating edge)
    {
        e0_t0 << -1.0, 0.0;
        e1_t0 << 1.0, 0.0;
        p_t0 << 0.0, 0.5;

        de0 << 1.6730970740318298, 0.8025388419628143;
        de1 << -1.616142749786377, -0.6420311331748962;
        dp << 0.0, -1.0;

        expected_toi = 0.4482900963;
    }

    SECTION("random")
    {
        // This times out on Windows debug
#if defined(NDEBUG) || !(defined(WIN32) || defined(_WIN32) || defined(__WIN32))
        auto impact = GENERATE(random_impacts(100, /*rigid=*/true));
#else
        auto impact = GENERATE(random_impacts(10, /*rigid=*/true));
#endif

        p_t0 = impact.p_t0;
        e0_t0 = impact.e0_t0;
        e1_t0 = impact.e1_t0;

        dp = impact.dp;
        de0 = impact.de0;
        de1 = impact.de1;

        expected_toi = impact.toi;
    }

#ifdef NDEBUG           // This case takes forever to run
    SECTION("parallel") //  (alpha=0 || alpha = 1)
    {
        p_t0 << 0.5, 0.0;
        e0_t0 << -0.5, 0.0;
        e1_t0 << -1.5, 0.0;

        // touches, intersects, passes-trough
        // double vel = GENERATE(1.0 + 1e-6, 2.0, 4.0);
        double vel = 2.0;
        // moving: both (same), ij, both (op), kl, both (same)
        // double j = GENERATE(0, 1, 2, 3, 4);
        double j = 0;
        // extension, no-deform, compression,
        // double dy = GENERATE(0.5, 0.0, -0.5);
        double dy = 0.0;

        dp << -(3 - j) * vel / 2.0, 0.0;
        de0 << (j - 1.0) * vel / 2.0, 0.0;
        // we only move one so we don't change the toi
        de1 << (j - 1.0) * vel / 2.0, dy;

        expected_toi = 1.0 / vel;
    }
#endif

    VectorMax3d p_t1 = p_t0 + dp;
    VectorMax3d e0_t1 = e0_t0 + de0;
    VectorMax3d e1_t1 = e1_t0 + de1;

    if (dim == 3) {
        p_t0.conservativeResize(3);
        p_t1.conservativeResize(3);
        e0_t0.conservativeResize(3);
        e0_t1.conservativeResize(3);
        e1_t0.conservativeResize(3);
        e1_t1.conservativeResize(3);
        p_t0[2] = p_t1[2] = e0_t0[2] = e1_t0[2] = e0_t1[2] = e1_t1[2] = 0;
    }

    check_toi(p_t0, e0_t0, e1_t0, p_t1, e0_t1, e1_t1, expected_toi);
    // Flip the edges
    check_toi(p_t0, e1_t0, e0_t0, p_t1, e1_t1, e0_t1, expected_toi);
}