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
        check_lanes(
            "coordinates", edge_edge_closest_point(a, b, c, d),
            [&](int l) {
                return edge_edge_closest_point(A[l], B[l], C[l], D[l]).eval();
            },
            VALUE_TOL);
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

#endif
