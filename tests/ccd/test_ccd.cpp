#include <catch2/catch.hpp>

#include <ipc/ipc.hpp>
#include <ipc/ccd/ccd.hpp>

#include <test_utils.hpp>

#include "collision_generator.hpp"

using namespace ipc;

static const double EPSILON = std::numeric_limits<float>::epsilon();

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

    Eigen::Vector2d p_t1 = p_t0 + dp;
    Eigen::Vector2d e0_t1 = e0_t0 + de0;
    Eigen::Vector2d e1_t1 = e1_t0 + de1;

    check_toi(p_t0, e0_t0, e1_t0, p_t1, e0_t1, e1_t1, expected_toi);
    // Flip the edges
    check_toi(p_t0, e1_t0, e0_t0, p_t1, e1_t1, e0_t1, expected_toi);
}

TEST_CASE("Repeated CCD", "[ccd][repeat]")
{
    const double FIRST_TOL = 1e-6, SECOND_TOL = 1e-7;
    const double FIRST_MAX_ITER = 1e6, SECOND_MAX_ITER = 1e6;

    // BroadPhaseMethod broadphase_method =
    //     GENERATE(BroadPhaseMethod::HASH_GRID, BroadPhaseMethod::BRUTE_FORCE);
    BroadPhaseMethod broadphase_method = BroadPhaseMethod::HASH_GRID;
    bool ignore_codimensional_vertices = true;
    double inflation_radius = 0;

    bool recompute_candidates = GENERATE(false, true);

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
    SECTION("hip small repeated toi")
    {
        t0_filename = "ccd-failure/small_repeated_toi_hip_0.obj";
        t1_filename = "ccd-failure/small_repeated_toi_hip_1.obj";
    }
    SECTION("hip inf-repeat 0")
    {
        t0_filename = "ccd-failure/inf_repeated_toi_hip_0.obj";
        t1_filename = "ccd-failure/inf_repeated_toi_hip_1.obj";
    }
    SECTION("hip inf-repeat 1")
    {
        t0_filename = "ccd-failure/s0121.obj";
        t1_filename = "ccd-failure/s1121.obj";
    }
    SECTION("hip inf-repeat 2")
    {
        t0_filename = "ccd-failure/s0110.obj";
        t1_filename = "ccd-failure/s1110.obj";
    }

    Eigen::MatrixXd V0, V1;
    Eigen::MatrixXi E, F;
    bool success =
        load_mesh(t0_filename, V0, E, F) && load_mesh(t1_filename, V1, E, F);
    if (!success) {
        return;
    }
    // REQUIRE(success);

    Candidates candidates;
    construct_collision_candidates(
        V0, V1, E, F, candidates, inflation_radius, broadphase_method,
        ignore_codimensional_vertices);

    bool has_collisions = !is_step_collision_free(
        candidates, V0, V1, E, F, FIRST_TOL, FIRST_MAX_ITER);

    double stepsize = compute_collision_free_stepsize(
        candidates, V0, V1, E, F, FIRST_TOL, FIRST_MAX_ITER);

    if (!has_collisions) {
        CHECK(stepsize == 1.0);
        return;
    }

    double collision_free_step_size = stepsize;
    bool has_collisions_repeated;
    double stepsize_repeated;
    do {
        Eigen::MatrixXd Vt = (V1 - V0) * collision_free_step_size + V0;
        // CHECK(!has_intersections(Vt, E, F));

        if (recompute_candidates) {
            construct_collision_candidates(
                V0, Vt, E, F, candidates, inflation_radius, broadphase_method,
                ignore_codimensional_vertices);
        }

        has_collisions_repeated = !is_step_collision_free(
            candidates, V0, Vt, E, F, SECOND_TOL, SECOND_MAX_ITER);

        stepsize_repeated = compute_collision_free_stepsize(
            candidates, V0, Vt, E, F, SECOND_TOL, SECOND_MAX_ITER);

        CAPTURE(
            t0_filename, t1_filename, broadphase_method, recompute_candidates,
            has_collisions, collision_free_step_size, has_collisions_repeated,
            stepsize_repeated);
        CHECK(!has_collisions_repeated);
        CHECK(stepsize_repeated == 1.0);

        collision_free_step_size *= stepsize_repeated;
    } while (has_collisions_repeated && stepsize_repeated != 1.0);
}

// This test takes too long
// TEST_CASE("Point-Triangle 3D CCD", "[ccd][point-triangle][timeout]")
// {
//     // point
//     double v0z = GENERATE(0.0, -1.0);
//     Eigen::Vector3d v0(0, 1, v0z);
//     // triangle = (v1, v2, v3)
//     Eigen::Vector3d v1(-1, 0, 1);
//     Eigen::Vector3d v2(1, 0, 1);
//     Eigen::Vector3d v3(0, 0, -1);
//
//     // displacements
//     double u0y =
//         -GENERATE(-1.0, 0.0, 0.5 - EPSILON, 0.5, 0.5 + EPSILON, 1.0, 2.0);
//     double u0z = GENERATE(-EPSILON, 0.0, EPSILON);
//     Eigen::Vector3d u0(0, u0y, u0z);
//     double u1y =
//         GENERATE(-1.0, 0.0, 0.5 - EPSILON, 0.5, 0.5 + EPSILON, 1.0, 2.0);
//     Eigen::Vector3d u1(0, u1y, 0);
//
//     bool is_collision_expected = ((-u0y + u1y >= 1) && (v0z + u0z >=
//     v3.z()));
//
//     double toi;
//     bool is_colliding = point_triangle_ccd(
//         v0, v1, v2, v3, v0 + u0, v1 + u1, v2 + u1, v3 + u1, toi);
//
//     CAPTURE(v0z, u0y, u1y, u0z, EPSILON);
//     REQUIRE(is_colliding >= is_collision_expected);
// }

// Disable these tests in debug because they take a long time.
#ifdef NDEBUG

TEST_CASE("Edge-Edge 3D CCD", "[ccd][edge-edge]")
{
    // e0 = (v0, v1)
    Eigen::Vector3d v0(-1, -1, 0);
    Eigen::Vector3d v1(1, -1, 0);
    // e2 = (v2, v3)
    // double e1x = GENERATE(
    //     -1 - EPSILON, -1, -1 + EPSILON, -0.5, 0, 0.5, 1 - EPSILON, 1,
    //     1 + EPSILON);
    double e1x = GENERATE(0);
    Eigen::Vector3d v2(e1x, 1, -1);
    Eigen::Vector3d v3(e1x, 1, 1);

    // displacements
    double y_displacement =
        GENERATE(-1.0, 0.0, 1 - EPSILON, 1.0, 1 + EPSILON, 2.0);

    Eigen::Vector3d u0, u1;
    bool is_collision_expected;
    SECTION("moving")
    {
        u0 << 0, y_displacement, 0;
        u1 << 0, -y_displacement, 0;
        is_collision_expected = y_displacement >= 1.0 && e1x >= -1 && e1x <= 1;
    }
    SECTION("fixed")
    {
        u0 << 0, 2 * y_displacement, 0;
        u1.setZero();
        is_collision_expected = y_displacement >= 2.0 && e1x >= -1 && e1x <= 1;
    }

    double toi;
    bool is_colliding =
        edge_edge_ccd(v0, v1, v2, v3, v0 + u0, v1 + u0, v2 + u1, v3 + u1, toi);

    CAPTURE(y_displacement, e1x);
    CHECK(is_colliding >= is_collision_expected);
}

// Only run this slow case on Linux and macOS
#if !defined(WIN32)

TEST_CASE("Zhongshi test case", "[ccd][point-triangle][zhongshi]")
{
    double qy = GENERATE(-EPSILON, 0, EPSILON);

    Eigen::Vector3d q;
    q << 0, qy, 0;

    Eigen::Vector3d b0;
    b0 << 0, 0, 0;
    Eigen::Vector3d b1;
    b1 << 0, 1, 0;
    Eigen::Vector3d b2;
    b2 << 1, 0, 0;

    Eigen::Vector3d t0;
    t0 << 0, 0, 1;
    Eigen::Vector3d t1;
    t1 << 0, 1, 1;
    Eigen::Vector3d t2;
    t2 << 1, 0, 1;

    Eigen::Vector3d q1;
    q1 << 0, qy, 0;

    bool is_collision_expected = q.y() >= 0;

    double toi;
    bool is_colliding = point_triangle_ccd(q, b0, b1, b2, q1, t0, t1, t2, toi);

    CAPTURE(qy);
    CHECK(is_colliding >= is_collision_expected);
}

#endif

TEST_CASE("Bolun test case", "[ccd][point-triangle][bolun]")
{
    Eigen::Vector3d x0(0.1, 0.1, 0.1), x1(0, 0, 1), x2(1, 0, 1), x3(0, 1, 1),
        x0b(0.1, 0.1, 0.1), x1b(0, 0, 0), x2b(0, 1, 0), x3b(1, 0, 0);

    bool is_collision_expected = true;

    double toi;
    bool is_colliding =
        point_triangle_ccd(x0, x1, x2, x3, x0b, x1b, x2b, x3b, toi);

    CHECK(is_colliding >= is_collision_expected);
}

TEST_CASE("Dobule root test case", "[ccd][edge-edge][double-root]")
{
    const Eigen::Vector3d a0s(-30022200, 2362580, 165247);
    const Eigen::Vector3d a1s(-32347850, 8312380, -1151003);
    const Eigen::Vector3d a0e(-28995600, 345838, 638580);
    const Eigen::Vector3d a1e(-31716930, 6104858, -713340);
    const Eigen::Vector3d b0(-30319900, 3148750, 0);
    const Eigen::Vector3d b1(-28548800, 900349, 0);

    bool is_collision_expected = true;

    double toi;
    bool is_colliding = edge_edge_ccd(a0s, a1s, b0, b1, a0e, a1e, b0, b1, toi);

    CHECK(is_colliding >= is_collision_expected);
}

TEST_CASE("Double root test case 2", "[ccd][edge-edge][double-root]")
{
    const Eigen::Vector3d a0s(0, 0, 1);
    const Eigen::Vector3d a1s(0, 1, 1);
    Eigen::Vector3d a0e(1, 1, 0);
    Eigen::Vector3d a1e(0, 0, 0);
    const Eigen::Vector3d b0(0.1, 0.2, 2);
    const Eigen::Vector3d b1(0.1, 0.2, -1);

    double t = GENERATE(0.5, 0.8, 0.88, 0.9, 1.0);
    a0e = (a0e - a0s) * t + a0s;
    a1e = (a1e - a1s) * t + a1s;

    bool is_collision_expected = true;

    double toi;
    bool is_colliding = edge_edge_ccd(a0s, a1s, b0, b1, a0e, a1e, b0, b1, toi);

    CHECK(is_colliding >= is_collision_expected);
}

#endif

TEST_CASE("Point-Plane CCD", "[ccd][point-plane]")
{
    Eigen::Vector3d p_t0 = Eigen::Vector3d::Random();

    Eigen::Vector3d origin = Eigen::Vector3d::Random();
    origin.x() = origin.z() = 0;
    Eigen::Vector3d normal(0, 1, 0);

    Eigen::Vector3d p_t1 = p_t0;
    double t = 0;
    SECTION("Known t")
    {
        t = GENERATE(1e-16, 1e-8, 1e-6, 1e-4, 0.1, 0.5, 0.8, 0.88, 0.9, 1.0);
    }
    SECTION("Random t") { t = GENERATE(take(100, random(0.0, 1.0))); }
    if (t == 0) {
        return;
    }
    p_t1.y() = (origin.y() + (t - 1) * p_t0.y()) / t;

    double toi;
    bool is_colliding =
        point_static_plane_ccd(p_t0, p_t1, origin, normal, toi, 0.99);

    CAPTURE(p_t0.x(), p_t0.y(), p_t0.z(), p_t1.x(), p_t1.y(), p_t1.z());
    CHECK(is_colliding);
    CHECK(toi <= t);
}
