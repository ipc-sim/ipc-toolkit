#include <catch2/catch.hpp>

#include <ipc/ipc.hpp>
#include <ipc/ccd/ccd.hpp>

#include <test_utils.hpp>

#include "collision_generator.hpp"

using namespace ipc;

TEST_CASE("Point-edge 2D CCD", "[ccd]")
{
    Eigen::Vector2d p_t0, p_t1, e0_t0, e0_t1, e1_t0, e1_t1;
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

    double toi, alpha;
    bool is_colliding = point_edge_ccd(
        p_t0, e0_t0, e1_t0, p_t1, e0_t1, e1_t1, toi,
        /*tmax=*/1.0,
        /*tolerance=*/1e-6,
        /*max_iterations=*/1e7,
        /*conservative_rescaling=*/1.0);
    REQUIRE(is_colliding == is_collision_expected);
    if (is_collision_expected) {
        CAPTURE(
            p_t0.transpose(), e0_t0.transpose(), e1_t0.transpose(),
            p_t1.transpose(), e0_t1.transpose(), e1_t1.transpose());
        CHECK(toi <= toi_expected);
    }
}

// ---------------------------------------------------
// Helper functions
// ---------------------------------------------------

/// Compares the time of impact of different implementations
/// against the expected time of impact
void check_toi(
    const Eigen::Vector2d& p_t0,
    const Eigen::Vector2d& e0_t0,
    const Eigen::Vector2d& e1_t0,
    const Eigen::Vector2d& p_t1,
    const Eigen::Vector2d& e0_t1,
    const Eigen::Vector2d& e1_t1,
    const double toi_expected)
{
    double toi_actual;
    // check autodiff code
    toi_actual = -1.0;
    bool has_collision = point_edge_ccd(
        p_t0, e0_t0, e1_t0, p_t1, e0_t1, e1_t1, toi_actual,
        /*tmax=*/1.0,
        /*tolerance=*/1e-6,
        /*max_iterations=*/1e7,
        /*conservative_rescaling=*/1.0);
    CHECK(has_collision);
    CHECK(toi_actual <= toi_expected);
}

// ---------------------------------------------------
// Tests
// ---------------------------------------------------
#if defined(NDEBUG) || !(defined(WIN32) || defined(_WIN32) || defined(__WIN32))

TEST_CASE("Point-edge 2D ToI", "[ccd][toi]")
{
    Eigen::Vector2d p_t0, e0_t0, e1_t0;
    Eigen::Vector2d dp, de0, de1;

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
        auto impact = GENERATE(random_impacts(100, /*rigid=*/true));

        p_t0 = impact.p_t0;
        e0_t0 = impact.e0_t0;
        e1_t0 = impact.e1_t0;

        dp = impact.dp;
        de0 = impact.de0;
        de1 = impact.de1;

        expected_toi = impact.toi;
    }

#ifdef NDEBUG
    SECTION("tangential") //  (alpha=0 || alpha = 1)
    {
        p_t0 << 0.5, 0.0;
        e0_t0 << -0.5, 0.0;
        e1_t0 << -1.5, 0.0;

        // touches, intersects, passes-trough
        double vel = GENERATE(1.0 + 1e-6, 2.0, 4.0);
        // moving: both (same), ij, both (op), kl, both (same)
        double j = GENERATE(0, 1, 2, 3, 4);
        // extension, no-deform, compression,
        double dy = GENERATE(0.5, 0.0, -0.5);

        dp << -(3 - j) * vel / 2.0, 0.0;
        de0 << (j - 1.0) * vel / 2.0, 0.0;
        // we only move one so we don't change the toi
        de1 << (j - 1.0) * vel / 2.0, dy;

        expected_toi = 1.0 / vel;
    }
#endif

    Eigen::Vector2d p_t1 = p_t0 + dp;
    Eigen::Vector2d e0_t1 = e0_t0 + de0;
    Eigen::Vector2d e1_t1 = e1_t0 + de1;

    check_toi(p_t0, e0_t0, e1_t0, p_t1, e0_t1, e1_t1, expected_toi);
    // Flip the edges
    check_toi(p_t0, e1_t0, e0_t0, p_t1, e1_t1, e0_t1, expected_toi);
}

#endif

TEST_CASE("Repeated CCD", "[ccd][thisone]")
{
    BroadPhaseMethod broadphase_method =
        GENERATE(BroadPhaseMethod::HASH_GRID, BroadPhaseMethod::BRUTE_FORCE);

    std::string t0_filename, t1_filename;
    SECTION("tooth")
    {
        t0_filename = "ccd-failure/repeated_toi_tooth_0.obj";
        t1_filename = "ccd-failure/repeated_toi_tooth_1.obj";
    }
    SECTION("hip")
    {
        t0_filename = "ccd-failure/repeated_toi_hip_0.obj";
        t1_filename = "ccd-failure/repeated_toi_hip_1.obj";
    }

    Eigen::MatrixXd V0, V1;
    Eigen::MatrixXi E, F;
    bool success =
        load_mesh(t0_filename, V0, E, F) && load_mesh(t1_filename, V1, E, F);
    if (!success) {
        return;
    }
    // REQUIRE(success);

    // REQUIRE(!has_intersections(V0, E, F));

    bool has_collisions =
        !is_step_collision_free(V0, V1, E, F, broadphase_method);

    if (!has_collisions) {
        return;
    }

    double stepsize =
        compute_collision_free_stepsize(V0, V1, E, F, broadphase_method);

    double collision_free_step_size = stepsize;
    bool has_collisions_repeated;
    double stepsize_repeated;
    do {
        Eigen::MatrixXd Vt = (V1 - V0) * collision_free_step_size + V0;
        CHECK(!has_intersections(Vt, E, F));

        has_collisions_repeated =
            !is_step_collision_free(V0, Vt, E, F, broadphase_method);

        stepsize_repeated =
            compute_collision_free_stepsize(V0, Vt, E, F, broadphase_method);

        CAPTURE(
            t0_filename, t1_filename, broadphase_method, has_collisions,
            collision_free_step_size, has_collisions_repeated,
            stepsize_repeated);
        CHECK(!has_collisions_repeated);
        CHECK(stepsize_repeated == 1.0);

        collision_free_step_size *= stepsize_repeated;
    } while (has_collisions_repeated);
}
