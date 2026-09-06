#include <catch2/catch_test_macros.hpp>

#include <tests/simd_utils.hpp>

#ifdef IPC_TOOLKIT_WITH_SIMD

#include <ipc/distance/edge_edge_mollifier.hpp>

#include <Eigen/Geometry>

#include <array>
#include <cmath>
#include <string>

using namespace ipc;
using namespace ipc::tests;

namespace {

/// @brief The threshold of two unit-length rest edges.
constexpr double EPS_X = 1e-3;

/// @brief Angles between the two edges placing a lane on either side of the
/// mollifier threshold: exactly parallel, mollified, just inside the
/// threshold, and inactive. The mollifier activates at sin²θ < eps_x, i.e.
/// θ < 0.0316 rad for unit rest edges.
constexpr std::array<double, 4> THETAS = { { 0.0, 0.02, 0.0316, 0.4 } };

/// @brief A lane-dependent rotation, so the gradients are not mostly zeros the
/// way an axis-aligned configuration would make them.
Eigen::Matrix3d rotation(const int l)
{
    return (Eigen::AngleAxisd(
                0.7 + 0.31 * l, Eigen::Vector3d(1, 2, 3).normalized())
            * Eigen::AngleAxisd(
                0.4 * l - 0.2, Eigen::Vector3d(-2, 1, 0.5).normalized()))
        .toRotationMatrix();
}

struct Config {
    Eigen::Vector3d ea0, ea1, eb0, eb1;
};

/// @brief Two unit-length edges meeting at an angle `theta`, rigidly placed by
/// lane. Unit lengths keep the rest threshold at exactly `EPS_X`, so `theta`
/// alone decides which side of the branch the lane lands on.
Config config(const double theta, const int l)
{
    const Eigen::Matrix3d R = rotation(l);
    const Eigen::Vector3d t(0.1 * l, -0.2, 0.3);

    const Eigen::Vector3d ea0(0, 0, 0);
    const Eigen::Vector3d ea1(1, 0, 0);
    const Eigen::Vector3d eb0(0, 1, 0.3);
    const Eigen::Vector3d eb1 =
        eb0 + Eigen::Vector3d(std::cos(theta), std::sin(theta), 0);

    return { R * ea0 + t, R * ea1 + t, R * eb0 + t, R * eb1 + t };
}

} // namespace

TEST_CASE(
    "SIMD batch edge-edge mollifier matches the scalar one lane-wise",
    "[distance][mollifier][simd]")
{
    // x straddling eps_x, so one batch carries mollified and inactive lanes.
    constexpr std::array<double, 4> REGION_XS = {
        { 0.0, 0.1 * EPS_X, 0.9 * EPS_X, 2.0 * EPS_X }
    };

    const auto check = [&](const std::string& name, auto&& scalar_fn,
                           auto&& batch_fn) {
        check_swept_lanes(
            name, REGION_XS, [&](double x) { return scalar_fn(x, EPS_X); },
            [&](Batch x) { return batch_fn(x, Batch(EPS_X)); });
    };

    check(
        "mollifier",
        [](double x, double e) { return edge_edge_mollifier(x, e); },
        [](Batch x, Batch e) { return edge_edge_mollifier(x, e); });
    check(
        "gradient",
        [](double x, double e) { return edge_edge_mollifier_gradient(x, e); },
        [](Batch x, Batch e) { return edge_edge_mollifier_gradient(x, e); });
    check(
        "hessian",
        [](double x, double e) { return edge_edge_mollifier_hessian(x, e); },
        [](Batch x, Batch e) { return edge_edge_mollifier_hessian(x, e); });
    check(
        "derivative_wrt_eps_x",
        [](double x, double e) {
            return edge_edge_mollifier_derivative_wrt_eps_x(x, e);
        },
        [](Batch x, Batch e) {
            return edge_edge_mollifier_derivative_wrt_eps_x(x, e);
        });
    check(
        "gradient_derivative_wrt_eps_x",
        [](double x, double e) {
            return edge_edge_mollifier_gradient_derivative_wrt_eps_x(x, e);
        },
        [](Batch x, Batch e) {
            return edge_edge_mollifier_gradient_derivative_wrt_eps_x(x, e);
        });
}

TEST_CASE(
    "SIMD batch edge-edge mollifier blends the threshold branch per lane",
    "[distance][mollifier][simd]")
{
    for (int offset = 0; offset < int(THETAS.size()); ++offset) {
        const Lanes thetas = lane_cases(THETAS, offset);
        const std::string at = " (offset " + std::to_string(offset) + ")";

        Points<3> EA0, EA1, EB0, EB1;
        // The rest configuration is the same geometry at theta = 0: both edges
        // are unit length, so the threshold is exactly EPS_X on every lane.
        Points<3> RA0, RA1, RB0, RB1;
        for (int l = 0; l < L; ++l) {
            const Config c = config(thetas[l], l);
            EA0[l] = c.ea0, EA1[l] = c.ea1, EB0[l] = c.eb0, EB1[l] = c.eb1;
            const Config r = config(0.0, l + 1);
            RA0[l] = r.ea0, RA1[l] = r.ea1, RB0[l] = r.eb0, RB1[l] = r.eb1;
        }

        const Eigen::Vector3<Batch> ea0 = pack(EA0), ea1 = pack(EA1),
                                    eb0 = pack(EB0), eb1 = pack(EB1);
        const Eigen::Vector3<Batch> ra0 = pack(RA0), ra1 = pack(RA1),
                                    rb0 = pack(RB0), rb1 = pack(RB1);
        const Batch eps_x(EPS_X);

        // The threshold itself must agree before anything built on it does.
        check_scalar_lanes(
            "threshold" + at, edge_edge_mollifier_threshold(ra0, ra1, rb0, rb1),
            [&](int l) {
                return edge_edge_mollifier_threshold(
                    RA0[l], RA1[l], RB0[l], RB1[l]);
            });
        check_scalar_lanes(
            "cross_squarednorm" + at,
            edge_edge_cross_squarednorm(ea0, ea1, eb0, eb1), [&](int l) {
                return edge_edge_cross_squarednorm(
                    EA0[l], EA1[l], EB0[l], EB1[l]);
            });
        check_scalar_lanes(
            "mollifier" + at, edge_edge_mollifier(ea0, ea1, eb0, eb1, eps_x),
            [&](int l) {
                return edge_edge_mollifier(
                    EA0[l], EA1[l], EB0[l], EB1[l], EPS_X);
            });

        check_lanes(
            "cross_squarednorm_gradient" + at,
            edge_edge_cross_squarednorm_gradient(ea0, ea1, eb0, eb1),
            [&](int l) {
                return edge_edge_cross_squarednorm_gradient(
                    EA0[l], EA1[l], EB0[l], EB1[l]);
            });
        check_lanes(
            "cross_squarednorm_hessian" + at,
            edge_edge_cross_squarednorm_hessian(ea0, ea1, eb0, eb1),
            [&](int l) {
                return edge_edge_cross_squarednorm_hessian(
                    EA0[l], EA1[l], EB0[l], EB1[l]);
            });

        check_lanes(
            "mollifier_gradient" + at,
            edge_edge_mollifier_gradient(ea0, ea1, eb0, eb1, eps_x),
            [&](int l) {
                return edge_edge_mollifier_gradient(
                    EA0[l], EA1[l], EB0[l], EB1[l], EPS_X);
            });
        check_lanes(
            "mollifier_hessian" + at,
            edge_edge_mollifier_hessian(ea0, ea1, eb0, eb1, eps_x), [&](int l) {
                return edge_edge_mollifier_hessian(
                    EA0[l], EA1[l], EB0[l], EB1[l], EPS_X);
            });

        check_lanes(
            "threshold_gradient" + at,
            edge_edge_mollifier_threshold_gradient(ra0, ra1, rb0, rb1),
            [&](int l) {
                return edge_edge_mollifier_threshold_gradient(
                    RA0[l], RA1[l], RB0[l], RB1[l]);
            });
        check_lanes(
            "mollifier_gradient_wrt_x" + at,
            edge_edge_mollifier_gradient_wrt_x(
                ra0, ra1, rb0, rb1, ea0, ea1, eb0, eb1),
            [&](int l) {
                return edge_edge_mollifier_gradient_wrt_x(
                    RA0[l], RA1[l], RB0[l], RB1[l], EA0[l], EA1[l], EB0[l],
                    EB1[l]);
            });
        check_lanes(
            "mollifier_gradient_jacobian_wrt_x" + at,
            edge_edge_mollifier_gradient_jacobian_wrt_x(
                ra0, ra1, rb0, rb1, ea0, ea1, eb0, eb1),
            [&](int l) {
                return edge_edge_mollifier_gradient_jacobian_wrt_x(
                    RA0[l], RA1[l], RB0[l], RB1[l], EA0[l], EA1[l], EB0[l],
                    EB1[l]);
            });
    }
}

#endif
