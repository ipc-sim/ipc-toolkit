#include <catch2/catch_test_macros.hpp>

#include <tests/simd_utils.hpp>

#ifdef IPC_TOOLKIT_WITH_SIMD

#include <ipc/tangent/relative_velocity.hpp>

using namespace ipc;
using namespace ipc::tests;

namespace {

// A relative velocity is a linear combination of the input velocities, so
// unlike a distance derivative it involves no cancellation of large terms and
// the two paths should agree to a rounding step. Hence VALUE_TOL rather than
// the DERIVATIVE_TOL that check_lanes defaults to.
constexpr double TOL = VALUE_TOL;

/// @brief A distinct barycentric coordinate per lane, so the Jacobians and
/// dx_dbeta tensors -- which depend on the coordinate and nothing else --
/// differ between lanes and the comparison is not vacuous.
Lanes lane_alphas()
{
    Lanes alphas;
    for (int l = 0; l < L; ++l) {
        alphas[l] = 0.15 + 0.2 * l;
    }
    return alphas;
}

Points<2> lane_coords()
{
    Points<2> coords;
    for (int l = 0; l < L; ++l) {
        coords[l] = Eigen::Vector2d(0.2 + 0.15 * l, 0.35 - 0.1 * l);
    }
    return coords;
}

} // namespace

TEST_CASE(
    "SIMD batch relative velocities match the scalar ones lane-wise (3D)",
    "[relative_velocity][simd]")
{
    const Points<3> A = random_points<3>(41), B = random_points<3>(42),
                    C = random_points<3>(43), D = random_points<3>(44);
    const Eigen::Vector3<Batch> a = pack(A), b = pack(B), c = pack(C),
                                d = pack(D);

    const Lanes ALPHA = lane_alphas();
    const Batch alpha = ALPHA.load();
    const Points<2> COORDS = lane_coords();
    const Eigen::Vector2<Batch> coords = pack(COORDS);

    SECTION("point-point")
    {
        check_lanes(
            "value", point_point_relative_velocity(a, b),
            [&](int l) {
                return point_point_relative_velocity(A[l], B[l]).eval();
            },
            TOL);
        check_lanes(
            "jacobian", point_point_relative_velocity_jacobian<Batch>(3),
            [&](int) {
                return point_point_relative_velocity_jacobian<double>(3);
            },
            TOL);
        // Γ does not depend on β here, so both paths must give exactly zero.
        check_lanes(
            "dx_dbeta", point_point_relative_velocity_dx_dbeta<Batch>(3),
            [&](int) {
                return point_point_relative_velocity_dx_dbeta<double>(3);
            },
            TOL);
    }

    SECTION("point-edge")
    {
        check_lanes(
            "value", point_edge_relative_velocity(a, b, c, alpha),
            [&](int l) {
                return point_edge_relative_velocity(A[l], B[l], C[l], ALPHA[l])
                    .eval();
            },
            TOL);
        check_lanes(
            "jacobian", point_edge_relative_velocity_jacobian(3, alpha),
            [&](int l) {
                return point_edge_relative_velocity_jacobian(3, ALPHA[l]);
            },
            TOL);
        check_lanes(
            "dx_dbeta", point_edge_relative_velocity_dx_dbeta(3, alpha),
            [&](int l) {
                return point_edge_relative_velocity_dx_dbeta(3, ALPHA[l]);
            },
            TOL);
    }

    SECTION("edge-edge")
    {
        check_lanes(
            "value", edge_edge_relative_velocity(a, b, c, d, coords),
            [&](int l) {
                return edge_edge_relative_velocity(
                           A[l], B[l], C[l], D[l], COORDS[l])
                    .eval();
            },
            TOL);
        check_lanes(
            "jacobian", edge_edge_relative_velocity_jacobian(coords),
            [&](int l) {
                return edge_edge_relative_velocity_jacobian(COORDS[l]).eval();
            },
            TOL);
        check_lanes(
            "dx_dbeta", edge_edge_relative_velocity_dx_dbeta(coords),
            [&](int l) {
                return edge_edge_relative_velocity_dx_dbeta(COORDS[l]).eval();
            },
            TOL);
    }

    SECTION("point-triangle")
    {
        check_lanes(
            "value", point_triangle_relative_velocity(a, b, c, d, coords),
            [&](int l) {
                return point_triangle_relative_velocity(
                           A[l], B[l], C[l], D[l], COORDS[l])
                    .eval();
            },
            TOL);
        check_lanes(
            "jacobian", point_triangle_relative_velocity_jacobian(coords),
            [&](int l) {
                return point_triangle_relative_velocity_jacobian(COORDS[l])
                    .eval();
            },
            TOL);
        check_lanes(
            "dx_dbeta", point_triangle_relative_velocity_dx_dbeta(coords),
            [&](int l) {
                return point_triangle_relative_velocity_dx_dbeta(COORDS[l])
                    .eval();
            },
            TOL);
    }
}

TEST_CASE(
    "SIMD batch relative velocities match the scalar ones lane-wise (2D)",
    "[relative_velocity][simd]")
{
    // Only point-point and point-edge are defined in 2D; the runtime `dim`
    // branch in the front ends is a size test, not a value test, so a batch
    // takes it the same way a scalar does.
    const Points<2> A = random_points<2>(51), B = random_points<2>(52),
                    C = random_points<2>(53);
    const Eigen::Vector2<Batch> a = pack(A), b = pack(B), c = pack(C);

    const Lanes ALPHA = lane_alphas();
    const Batch alpha = ALPHA.load();

    check_lanes(
        "point_point_value", point_point_relative_velocity(a, b),
        [&](int l) { return point_point_relative_velocity(A[l], B[l]).eval(); },
        TOL);
    check_lanes(
        "point_point_jacobian",
        point_point_relative_velocity_jacobian<Batch>(2),
        [&](int) { return point_point_relative_velocity_jacobian<double>(2); },
        TOL);
    check_lanes(
        "point_point_dx_dbeta",
        point_point_relative_velocity_dx_dbeta<Batch>(2),
        [&](int) { return point_point_relative_velocity_dx_dbeta<double>(2); },
        TOL);

    check_lanes(
        "point_edge_value", point_edge_relative_velocity(a, b, c, alpha),
        [&](int l) {
            return point_edge_relative_velocity(A[l], B[l], C[l], ALPHA[l])
                .eval();
        },
        TOL);
    check_lanes(
        "point_edge_jacobian", point_edge_relative_velocity_jacobian(2, alpha),
        [&](int l) {
            return point_edge_relative_velocity_jacobian(2, ALPHA[l]);
        },
        TOL);
    check_lanes(
        "point_edge_dx_dbeta", point_edge_relative_velocity_dx_dbeta(2, alpha),
        [&](int l) {
            return point_edge_relative_velocity_dx_dbeta(2, ALPHA[l]);
        },
        TOL);
}

#endif
