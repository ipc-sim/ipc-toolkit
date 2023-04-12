#include <catch2/catch_all.hpp>

#include <ipc/ipc.hpp>
#include <ipc/ccd/ccd.hpp>
#include <ipc/ccd/point_static_plane.hpp>

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
        /*min_distance=*/0.0,
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
        /*min_distance=*/0.0,
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
    const double MIN_DISTANCE = 0.0;

    // BroadPhaseMethod method = GENERATE_BROAD_PHASE_METHODS();
    BroadPhaseMethod broadphase_method = BroadPhaseMethod::HASH_GRID;
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

    CollisionMesh mesh = CollisionMesh::build_from_full_mesh(V0, E, F);
    // Discard codimensional/internal vertices
    V0 = mesh.vertices(V0);
    V1 = mesh.vertices(V1);

    Candidates candidates;
    candidates.build(mesh, V0, V1, inflation_radius, broadphase_method);

    bool has_collisions = !candidates.is_step_collision_free(
        mesh, V0, V1, MIN_DISTANCE, FIRST_TOL, FIRST_MAX_ITER);

    double stepsize = candidates.compute_collision_free_stepsize(
        mesh, V0, V1, MIN_DISTANCE, FIRST_TOL, FIRST_MAX_ITER);

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
            candidates.build(mesh, V0, Vt, inflation_radius, broadphase_method);
        }

        has_collisions_repeated = !candidates.is_step_collision_free(
            mesh, V0, Vt, MIN_DISTANCE, SECOND_TOL, SECOND_MAX_ITER);

        stepsize_repeated = candidates.compute_collision_free_stepsize(
            mesh, V0, Vt, MIN_DISTANCE, SECOND_TOL, SECOND_MAX_ITER);

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
    const Eigen::Vector3d a0s(-3.0022200, 0.2362580, 0.0165247);
    const Eigen::Vector3d a1s(-3.2347850, 0.8312380, -0.1151003);
    const Eigen::Vector3d a0e(-2.8995600, 0.0345838, 0.0638580);
    const Eigen::Vector3d a1e(-3.1716930, 0.6104858, -0.0713340);
    const Eigen::Vector3d b0(-3.0319900, 0.3148750, 0.0000000);
    const Eigen::Vector3d b1(-2.8548800, 0.0900349, 0.0000000);

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

TEST_CASE("No Zero ToI CCD", "[ccd][no-zero-toi]")
{
    Eigen::Vector3d p_t0, t0_t0, t1_t0, t2_t0, p_t1, t0_t1, t1_t1, t2_t1;
    p_t0 << 0.0133653, 0.100651, -0.0215935;
    t0_t0 << 0.0100485, 0.0950896, -0.0171013;
    t1_t0 << 0.0130388, 0.100666, -0.0218112;
    t2_t0 << 0.015413, 0.100554, -0.0202265;
    p_t1 << 0.0133652999767858, 0.099670000268615, -0.0215934999996444;
    t0_t1 << 0.0100484999799995, 0.0941086002577558, -0.0171012999972189;
    t1_t1 << 0.0130387999724314, 0.0996850002629403, -0.0218111999936902;
    t2_t1 << 0.0154129999740718, 0.0995730002646605, -0.020226499996014;

    double toi;
    bool is_impacting = point_triangle_ccd(
        p_t0, t0_t0, t1_t0, t2_t0, p_t1, t0_t1, t1_t1, t2_t1, toi);

    CAPTURE(toi);
    CHECK(!is_impacting);
}

TEST_CASE("Slow EE CCD", "[ccd][edge-edge][slow]")
{
    Eigen::Vector3d ea0_t0, ea1_t0, eb0_t0, eb1_t0, ea0_t1, ea1_t1, eb0_t1,
        eb1_t1;
    ea0_t0 << 1, 0.50803125, 2.10835646075301e-18;
    ea1_t0 << -2.38233935445388e-18, 0.50803125, 1;
    eb0_t0 << -4.99999999958867e-07, 0.5, 0;
    eb1_t0 << -4.99999999958867e-07, 0.5, 1;
    ea0_t1 << 1, 0.47124375, 4.11078309465837e-18;
    ea1_t1 << -2.8526707189104e-18, 0.47124375, 1;
    eb0_t1 << -4.99999999958867e-07, 0.5, 0;
    eb1_t1 << -4.99999999958867e-07, 0.5, 1;

    // BENCHMARK("compute toi")
    // {
    //     double toi;
    //     bool is_impacting = edge_edge_ccd(
    //         ea0_t0, ea1_t0, eb0_t0, eb1_t0, ea0_t1, ea1_t1, eb0_t1, eb1_t1,
    //         toi);
    // };

    double tol = DEFAULT_CCD_TOLERANCE;
    long max_iter = DEFAULT_CCD_MAX_ITERATIONS;
    while (tol < 1) {
        double toi;
        bool is_impacting = edge_edge_ccd(
            ea0_t0, ea1_t0, eb0_t0, eb1_t0, ea0_t1, ea1_t1, eb0_t1, eb1_t1, toi,
            0.0, 1.0, tol, max_iter);
        tol *= 10;

        CAPTURE(toi);
        CHECK(is_impacting);
    }
}

TEST_CASE("Slow EE CCD 2", "[ccd][edge-edge][slow][thisone]")
{
    Eigen::Vector3d ea0_t0, ea1_t0, eb0_t0, eb1_t0, ea0_t1, ea1_t1, eb0_t1,
        eb1_t1;
    ea0_t0 << 1.00002232466453, 0.500004786049044, -2.06727783590977e-05;
    ea1_t0 << 1.64687846177844e-05, 0.499996645067319, 1.63939999009028e-05;
    eb0_t0 << 1, 0.5, 0;
    eb1_t0 << 0, 0.5, 0;
    ea0_t1 << 1.00294282700155, 0.498652627047143, 0.003626320742036;
    ea1_t1 << -0.00219276550735626, 0.500871179186644, -0.00315828804921928;
    eb0_t1 << 1, 0.5, 0;
    eb1_t1 << 0, 0.5, 0;

    double tol = 1e-6;
    long max_iter = 1e6;
    double tmax = 2.8076171875e-03;

    double toi;
    bool is_impacting = edge_edge_ccd(
        ea0_t0, ea1_t0, eb0_t0, eb1_t0, ea0_t1, ea1_t1, eb0_t1, eb1_t1, toi,
        0.0, tmax, tol, max_iter);

    CAPTURE(toi);
    CHECK(is_impacting);
}

TEST_CASE("Squash Tet", "[ccd]")
{
    const double dhat = 1e-3;

    Eigen::MatrixXd rest_vertices(4, 3);
    rest_vertices << 0.0, 0.0, 0.0, //
        1.0, 0.0, 0.0,              //
        0.0, 1.0, 0.0,              //
        0.0, 0.0, 1.0;

    Eigen::MatrixXd deformed_vertices = rest_vertices;
    deformed_vertices(3, 0) += 0.1;
    deformed_vertices(3, 1) -= 0.1;
    deformed_vertices(3, 2) = -0.5 * dhat;

    Eigen::MatrixXi edges(6, 2);
    edges << 0, 1, //
        0, 2,      //
        0, 3,      //
        1, 2,      //
        1, 3,      //
        2, 3;
    Eigen::MatrixXi faces(4, 3);
    faces << 0, 2, 1, //
        0, 1, 3,      //
        0, 3, 2,      //
        1, 2, 3;

    ipc::CollisionMesh mesh =
        ipc::CollisionMesh::build_from_full_mesh(rest_vertices, edges, faces);

    // BENCHMARK("compute toi")
    // {
    logger().debug(
        "toi={}",
        ipc::compute_collision_free_stepsize(
            mesh, rest_vertices, deformed_vertices));
    // };
}