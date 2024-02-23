#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <ipc/config.hpp>
#include <ipc/ccd/ccd.hpp>
#include <ipc/ccd/additive_ccd.hpp>

using namespace ipc;

static const double EPSILON = std::numeric_limits<float>::epsilon();

#ifdef IPC_TOOLKIT_WITH_INEXACT_CCD
TEST_CASE("Edge-Edge CCD", "[ccd][3D][edge-edge][!mayfail]")
#else
TEST_CASE("Edge-Edge CCD", "[ccd][3D][edge-edge]")
#endif
{
    Eigen::Vector3d ea0_t0, ea1_t0, eb0_t0, eb1_t0, ea0_t1, ea1_t1, eb0_t1,
        eb1_t1;
    bool is_collision_expected;
    bool conservative_check = true;

    double tol = DEFAULT_CCD_TOLERANCE;
    long max_iter = DEFAULT_CCD_MAX_ITERATIONS;
    double tmax = 1;

#if !defined(WIN32) || defined(NDEBUG)
    SECTION("General")
    {
        double uy = GENERATE(-1.0, 0.0, 1 - EPSILON, 1.0, 1 + EPSILON, 2.0);
        double e1x = GENERATE(
            -1 - EPSILON, -1, -1 + EPSILON, -0.5, 0, 0.5, 1 - EPSILON, 1,
            1 + EPSILON);

        ea0_t0 << -1, -1, 0;
        ea1_t0 << 1, -1, 0;
        eb0_t0 << e1x, 1, -1;
        eb1_t0 << e1x, 1, 1;

        Eigen::Vector3d u0, u1;
        SECTION("moving")
        {
            u0 << 0, uy, 0;
            u1 << 0, -uy, 0;
            is_collision_expected = uy >= 1.0 && e1x >= -1 && e1x <= 1;
        }
        SECTION("fixed")
        {
            u0 << 0, 2 * uy, 0;
            u1.setZero();
            is_collision_expected = uy >= 2.0 && e1x >= -1 && e1x <= 1;
        }

        ea0_t1 = ea0_t0 + u0;
        ea1_t1 = ea1_t0 + u0;
        eb0_t1 = eb0_t0 + u1;
        eb1_t1 = eb1_t0 + u1;
    }
    SECTION("Double root test case 1")
    {
        ea0_t0 << -3.0022200, 0.2362580, 0.0165247;
        ea1_t0 << -3.2347850, 0.8312380, -0.1151003;
        eb0_t0 << -3.0319900, 0.3148750, 0.0000000;
        eb1_t0 << -2.8548800, 0.0900349, 0.0000000;
        ea0_t1 << -2.8995600, 0.0345838, 0.0638580;
        ea1_t1 << -3.1716930, 0.6104858, -0.0713340;
        eb0_t1 = eb0_t0;
        eb1_t1 = eb1_t0;

        is_collision_expected = true;
    }
    SECTION("Double root test case 2")
    {
        ea0_t0 << 0, 0, 1;
        ea1_t0 << 0, 1, 1;
        eb0_t0 << 0.1, 0.2, 2;
        eb1_t0 << 0.1, 0.2, -1;
        ea0_t1 << 1, 1, 0;
        ea1_t1 << 0, 0, 0;
        eb0_t1 = eb0_t0;
        eb1_t1 = eb1_t0;

        const double t = GENERATE(0.5, 0.8, 0.88, 0.9, 1.0);
        ea0_t1 = (ea0_t1 - ea0_t0) * t + ea0_t0;
        ea1_t1 = (ea1_t1 - ea1_t0) * t + ea1_t0;

        is_collision_expected = true;
    }
#endif
    SECTION("Slow Case 1")
    {
        ea0_t0 << 1, 0.50803125, 2.10835646075301e-18;
        ea1_t0 << -2.38233935445388e-18, 0.50803125, 1;
        eb0_t0 << -4.99999999958867e-07, 0.5, 0;
        eb1_t0 << -4.99999999958867e-07, 0.5, 1;
        ea0_t1 << 1, 0.47124375, 4.11078309465837e-18;
        ea1_t1 << -2.8526707189104e-18, 0.47124375, 1;
        eb0_t1 << -4.99999999958867e-07, 0.5, 0;
        eb1_t1 << -4.99999999958867e-07, 0.5, 1;

        tol = GENERATE(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1);

        is_collision_expected = true;
        conservative_check = false;
    }
    SECTION("Slow Case 2")
    {
        ea0_t0 << 1.00002232466453, 0.500004786049044, -2.06727783590977e-05;
        ea1_t0 << 1.64687846177844e-05, 0.499996645067319, 1.63939999009028e-05;
        eb0_t0 << 1, 0.5, 0;
        eb1_t0 << 0, 0.5, 0;
        ea0_t1 << 1.00294282700155, 0.498652627047143, 0.003626320742036;
        ea1_t1 << -0.00219276550735626, 0.500871179186644, -0.00315828804921928;
        eb0_t1 << 1, 0.5, 0;
        eb1_t1 << 0, 0.5, 0;

        max_iter = 1e6;
        tmax = 2.8076171875e-03;

        is_collision_expected = true;
        conservative_check = false;
    }
    SECTION("Adversarial ACCD Case")
    {
        const double d0 = GENERATE(1e-2, 1e-4, 1e-6, 1e-8);
        ea0_t0 << -1.0, -d0 / 2, 0.0;
        ea1_t0 << 1.0, -d0 / 2, 0.0;
        eb0_t0 << 0.0, d0 / 2, -1.0;
        eb1_t0 << 0.0, d0 / 2, 1.0;

        const double dy = d0;
        const double scale = 1e-3;

        ea0_t1 = ea0_t0 + Eigen::Vector3d(-scale, dy, -scale);
        ea1_t1 = ea1_t0 + Eigen::Vector3d(scale, dy, scale);
        eb0_t1 = eb0_t0 + Eigen::Vector3d(scale, -dy, scale);
        eb1_t1 = eb1_t0 + Eigen::Vector3d(-scale, -dy, -scale);

        // this ternary operator is to force MSVC to use 1 or 0
        is_collision_expected = dy >= d0 / 2;
        conservative_check = false;
    }
    CAPTURE(is_collision_expected);

    double toi;
    bool is_colliding = edge_edge_ccd(
        ea0_t0, ea1_t0, eb0_t0, eb1_t0, ea0_t1, ea1_t1, eb0_t1, eb1_t1, toi,
        /*min_distance=*/0.0, tmax, tol, max_iter);
    if (conservative_check) {
        CHECK((is_colliding || !is_collision_expected));
    } else {
        CHECK(is_colliding == is_collision_expected);
    }

    is_colliding = additive_ccd::edge_edge_ccd(
        ea0_t0, ea1_t0, eb0_t0, eb1_t0, ea0_t1, ea1_t1, eb0_t1, eb1_t1, toi,
        /*min_distance=*/0.0, tmax);
    if (conservative_check) {
        CHECK((is_colliding || !is_collision_expected));
    } else {
        CHECK(is_colliding == is_collision_expected);
    }
}
