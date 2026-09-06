#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_range.hpp>

#include <tests/simd_utils.hpp>

#ifdef IPC_TOOLKIT_WITH_SIMD

#include <ipc/tangent/closest_point.hpp>

#include <array>
#include <cmath>
#include <string>

using namespace ipc;
using namespace ipc::tests;

namespace {

/// @brief Two edges crossing at `theta`, offset along z so the closest points
/// land in the interior of both. That interior-interior case is the only one
/// these functions are called for, so it is what we test.
struct EdgePair {
    Eigen::Vector3d ea0, ea1, eb0, eb1;
};

EdgePair crossing_edges(const double theta)
{
    return { Eigen::Vector3d(-1, 0, 0), Eigen::Vector3d(1, 0, 0),
             Eigen::Vector3d(-std::cos(theta), -std::sin(theta), -0.5),
             Eigen::Vector3d(std::cos(theta), std::sin(theta), -0.5) };
}

/// @brief Two exactly parallel edges, separated in z and slid past each other
/// along their shared direction by `shift`.
///
/// Parallel is what makes the 2x2 Gram matrix singular. The shift is what
/// makes the case worth testing: without it the offset between the edges is
/// perpendicular to both, the right-hand side comes out exactly zero, and a
/// lane that was wrongly solved rather than zeroed would still leave no
/// residual to catch. Sliding the edges gives that lane a residual of order
/// ‖b‖, so the singular path has to actually be taken.
EdgePair parallel_edges(const double shift)
{
    return { Eigen::Vector3d(-1, 0, 0), Eigen::Vector3d(1, 0, 0),
             Eigen::Vector3d(-1 + shift, 0, -0.5),
             Eigen::Vector3d(1 + shift, 0, -0.5) };
}

} // namespace

TEST_CASE(
    "A batch of closest points agrees with the scalar path, one problem per "
    "lane",
    "[closest_point][simd]")
{
    const Points<3> A = random_points<3>(61), B = random_points<3>(62),
                    C = random_points<3>(63), D = random_points<3>(64);
    const Eigen::Vector3<Batch> a = pack(A), b = pack(B), c = pack(C),
                                d = pack(D);

    SECTION("point-edge")
    {
        check_scalar_lanes(
            "coordinate", point_edge_closest_point(a, b, c),
            [&](int l) { return point_edge_closest_point(A[l], B[l], C[l]); });
    }

    SECTION("point-triangle")
    {
        check_lanes(
            "coordinates", point_triangle_closest_point(a, b, c, d),
            [&](int l) {
                return point_triangle_closest_point(A[l], B[l], C[l], D[l])
                    .eval();
            },
            VALUE_TOL);
        check_lanes(
            "jacobian", point_triangle_closest_point_jacobian(a, b, c, d),
            [&](int l) {
                return point_triangle_closest_point_jacobian(
                    A[l], B[l], C[l], D[l]);
            });
    }

    SECTION("edge-edge")
    {
        // These coordinates get the derivative bound, not VALUE_TOL. The two
        // segments here are unconstrained in both direction and offset, so
        // the closest points routinely land far off the ends of both: a
        // coordinate of 20+ is an ordinary draw, not a degenerate one. Such a
        // coordinate is an amplification of the inputs, and the batch and
        // scalar paths contract multiply-adds differently, so their agreement
        // is bounded by the conditioning of the 2x2 solve rather than by a
        // rounding step on the answer. That is the same reason the
        // deliberately near-parallel test below uses this bound. A lane that
        // is structurally wrong still differs by order of its own magnitude,
        // far above this.
        check_lanes(
            "coordinates", edge_edge_closest_point(a, b, c, d),
            [&](int l) {
                return edge_edge_closest_point(A[l], B[l], C[l], D[l]).eval();
            },
            DERIVATIVE_TOL);
        check_lanes(
            "jacobian", edge_edge_closest_point_jacobian(a, b, c, d),
            [&](int l) {
                return edge_edge_closest_point_jacobian(A[l], B[l], C[l], D[l]);
            });
    }
}

TEST_CASE(
    "Closest points agree with the scalar path even for nearly parallel edges",
    "[closest_point][simd]")
{
    // The angle between the edges sets how well conditioned the 2x2 system is,
    // so we pack one batch with lanes ranging from well separated to nearly
    // parallel. That is the regime where the solve is worst conditioned, and
    // so where a batch and a scalar are most likely to drift apart. Rotating
    // the offset puts a different angle in each lane on every pass.
    constexpr std::array<double, 4> THETAS = { { 1.0, 0.1, 1e-2, 1e-3 } };
    const int offset = GENERATE(range(0, 4));

    const Lanes thetas = lane_cases(THETAS, offset);

    Points<3> EA0, EA1, EB0, EB1;
    for (int l = 0; l < L; ++l) {
        const EdgePair e = crossing_edges(thetas[l]);
        EA0[l] = e.ea0, EA1[l] = e.ea1, EB0[l] = e.eb0, EB1[l] = e.eb1;
    }

    check_lanes(
        "coordinates",
        edge_edge_closest_point(pack(EA0), pack(EA1), pack(EB0), pack(EB1)),
        [&](int l) {
            return edge_edge_closest_point(EA0[l], EA1[l], EB0[l], EB1[l])
                .eval();
        });
}

TEST_CASE(
    "A singular lane is zeroed without disturbing the lanes beside it",
    "[closest_point][simd]")
{
    // theta = 0 stands for a parallel pair, whose Gram matrix is singular, so
    // the solve zeroes that lane instead of dividing by the determinant. Every
    // offset below mixes such a lane with a well conditioned one, which is the
    // case a single batch has to answer two ways at once: zero one lane, solve
    // the other, and hold the residual check to the solved lane only.
    constexpr std::array<double, 4> THETAS = { { 0.0, 1.0, 0.0, 0.1 } };
    const int offset = GENERATE(range(0, 4));

    const Lanes thetas = lane_cases(THETAS, offset);

    Points<3> EA0, EA1, EB0, EB1;
    for (int l = 0; l < L; ++l) {
        const EdgePair e =
            thetas[l] == 0.0 ? parallel_edges(0.75) : crossing_edges(thetas[l]);
        EA0[l] = e.ea0, EA1[l] = e.ea1, EB0[l] = e.eb0, EB1[l] = e.eb1;
    }

    const Eigen::Vector2<Batch> coords =
        edge_edge_closest_point(pack(EA0), pack(EA1), pack(EB0), pack(EB1));

    check_lanes("coordinates", coords, [&](int l) {
        return edge_edge_closest_point(EA0[l], EA1[l], EB0[l], EB1[l]).eval();
    });

    // Pin down which lanes were the degenerate ones. Without this the check
    // above would still pass if the batch and the scalar path agreed on some
    // other answer for a parallel pair.
    for (int l = 0; l < L; ++l) {
        CAPTURE(l, thetas[l]);
        if (thetas[l] == 0.0) {
            CHECK(coords[0].get(l) == 0.0);
            CHECK(coords[1].get(l) == 0.0);
        } else {
            CHECK(coords[0].get(l) != 0.0);
        }
    }
}

#endif
