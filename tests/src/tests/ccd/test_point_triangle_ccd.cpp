#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <ipc/config.hpp>
#include <ipc/ccd/ccd.hpp>
#include <ipc/ccd/additive_ccd.hpp>

using namespace ipc;

static const double EPSILON = std::numeric_limits<float>::epsilon();

#ifdef IPC_TOOLKIT_WITH_INEXACT_CCD
TEST_CASE("Point-Triangle CCD", "[ccd][3D][point-triangle][!mayfail]")
#else
TEST_CASE("Point-Triangle CCD", "[ccd][3D][point-triangle]")
#endif
{
    Eigen::Vector3d p_t0, t0_t0, t1_t0, t2_t0, p_t1, t0_t1, t1_t1, t2_t1;
    bool is_collision_expected;
    bool conservative_check = true;

    std::string name;
    SECTION("General")
    {
        name = "General";
        const double v0z = GENERATE(0.0, -1.0);
        const double u0y =
            -GENERATE(-1.0, 0.0, 0.5 - EPSILON, 0.5, 0.5 + EPSILON, 1.0, 2.0);
        const double u0z = GENERATE(-EPSILON, 0.0, EPSILON);
        const double u1y =
            GENERATE(-1.0, 0.0, 0.5 - EPSILON, 0.5, 0.5 + EPSILON, 1.0, 2.0);
        CAPTURE(v0z, u0y, u1y, u0z);

        p_t0 << 0, 1, v0z;
        t0_t0 << -1, 0, 1;
        t1_t0 << 1, 0, 1;
        t2_t0 << 0, 0, -1;

        p_t1 = p_t0 + Eigen::Vector3d(0, u0y, u0z);
        t0_t1 = t0_t0 + Eigen::Vector3d(0, u1y, 0);
        t1_t1 = t1_t0 + Eigen::Vector3d(0, u1y, 0);
        t2_t1 = t2_t0 + Eigen::Vector3d(0, u1y, 0);

        is_collision_expected = ((-u0y + u1y >= 1) && (v0z + u0z >= t2_t0.z()));
    }
#ifdef NDEBUG
#ifndef WIN32
    SECTION("Zhongshi's test case")
    {
        name = "Zhongshi's test case";
        const double qy = GENERATE(-EPSILON, 0, EPSILON);
        CAPTURE(qy);

        p_t0 << 0, qy, 0;
        t0_t0 << 0, 0, 0;
        t1_t0 << 0, 1, 0;
        t2_t0 << 1, 0, 0;
        p_t1 << 0, qy, 0;
        t0_t1 << 0, 0, 1;
        t1_t1 << 0, 1, 1;
        t2_t1 << 1, 0, 1;

        is_collision_expected = qy >= 0;
    }
#endif
    SECTION("Bolun's test case")
    {
        name = "Bolun's test case";
        p_t0 << 0.1, 0.1, 0.1;
        t0_t0 << 0, 0, 1;
        t1_t0 << 1, 0, 1;
        t2_t0 << 0, 1, 1;
        p_t1 << 0.1, 0.1, 0.1;
        t0_t1 << 0, 0, 0;
        t1_t1 << 0, 1, 0;
        t2_t1 << 1, 0, 0;
        is_collision_expected = true;
    }
#endif
    SECTION("No Zero ToI CCD")
    {
        name = "No Zero ToI CCD";
        p_t0 << 0.0133653, 0.100651, -0.0215935;
        t0_t0 << 0.0100485, 0.0950896, -0.0171013;
        t1_t0 << 0.0130388, 0.100666, -0.0218112;
        t2_t0 << 0.015413, 0.100554, -0.0202265;
        p_t1 << 0.0133652999767858, 0.099670000268615, -0.0215934999996444;
        t0_t1 << 0.0100484999799995, 0.0941086002577558, -0.0171012999972189;
        t1_t1 << 0.0130387999724314, 0.0996850002629403, -0.0218111999936902;
        t2_t1 << 0.0154129999740718, 0.0995730002646605, -0.020226499996014;
        is_collision_expected = false;
        conservative_check = false;
    }
    INFO(name);

    double toi;
    bool is_colliding = point_triangle_ccd(
        p_t0, t0_t0, t1_t0, t2_t0, p_t1, t0_t1, t1_t1, t2_t1, toi);

    if (conservative_check) {
        CHECK((is_colliding || !is_collision_expected));
    } else {
        CHECK(is_colliding == is_collision_expected);
    }

    is_colliding = additive_ccd::point_triangle_ccd(
        p_t0, t0_t0, t1_t0, t2_t0, p_t1, t0_t1, t1_t1, t2_t1, toi);
    if (conservative_check) {
        CHECK((is_colliding || !is_collision_expected));
    } else {
        CHECK(is_colliding == is_collision_expected);
    }
}
