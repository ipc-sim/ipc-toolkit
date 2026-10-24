#include <catch2/catch_test_macros.hpp>

#include <tests/simd_utils.hpp>

#ifdef IPC_TOOLKIT_WITH_SIMD

#include <ipc/distance/signed/line_line.hpp>
#include <ipc/distance/signed/point_line.hpp>
#include <ipc/distance/signed/point_plane.hpp>

using namespace ipc;
using namespace ipc::tests;

namespace {

/// @brief Agreement bound for a signed distance *value*.
///
/// Looser than the shared `VALUE_TOL`, and structurally so: an unsigned
/// distance is a squared or absolute length, but a signed distance is
/// `normal.dot(p - x0)` with the normal a *normalized* cross product. The
/// cross product cancels, and the `sqrt` and division that follow amplify the
/// batch and scalar paths' differing multiply-add contraction rather than
/// costing one rounding step. Measured agreement is ~1e-10 relative.
constexpr double SIGNED_VALUE_TOL = DERIVATIVE_TOL;

} // namespace

TEST_CASE(
    "SIMD batch signed distances match the scalar ones lane-wise",
    "[distance][signed_distance][simd]")
{
    SECTION("point-line (2D)")
    {
        const auto P = random_points<2>(11);
        const auto E0 = random_points<2>(12);
        const auto E1 = random_points<2>(13);

        const auto p = pack(P), e0 = pack(E0), e1 = pack(E1);

        check_scalar_lanes(
            "value", point_line_signed_distance(p, e0, e1),
            [&](int l) {
                return point_line_signed_distance(P[l], E0[l], E1[l]);
            },
            SIGNED_VALUE_TOL);
        check_lanes(
            "gradient", point_line_signed_distance_gradient(p, e0, e1),
            [&](int l) {
                return point_line_signed_distance_gradient(P[l], E0[l], E1[l])
                    .eval();
            });
        check_lanes(
            "hessian", point_line_signed_distance_hessian(p, e0, e1),
            [&](int l) {
                return point_line_signed_distance_hessian(P[l], E0[l], E1[l])
                    .eval();
            });
    }

    SECTION("line-line")
    {
        const auto A = random_points<3>(21);
        const auto B = random_points<3>(22);
        const auto C = random_points<3>(23);
        const auto D = random_points<3>(24);

        const auto a = pack(A), b = pack(B), c = pack(C), d = pack(D);

        check_scalar_lanes(
            "value", line_line_signed_distance(a, b, c, d),
            [&](int l) {
                return line_line_signed_distance(A[l], B[l], C[l], D[l]);
            },
            SIGNED_VALUE_TOL);
        check_lanes(
            "gradient", line_line_signed_distance_gradient(a, b, c, d),
            [&](int l) {
                return line_line_signed_distance_gradient(
                           A[l], B[l], C[l], D[l])
                    .eval();
            });
        check_lanes(
            "hessian", line_line_signed_distance_hessian(a, b, c, d),
            [&](int l) {
                return line_line_signed_distance_hessian(A[l], B[l], C[l], D[l])
                    .eval();
            });
    }

    SECTION("point-plane")
    {
        const auto P = random_points<3>(31);
        const auto T0 = random_points<3>(32);
        const auto T1 = random_points<3>(33);
        const auto T2 = random_points<3>(34);

        const auto p = pack(P), t0 = pack(T0), t1 = pack(T1), t2 = pack(T2);

        check_scalar_lanes(
            "value", point_plane_signed_distance(p, t0, t1, t2),
            [&](int l) {
                return point_plane_signed_distance(P[l], T0[l], T1[l], T2[l]);
            },
            SIGNED_VALUE_TOL);
        check_lanes(
            "gradient", point_plane_signed_distance_gradient(p, t0, t1, t2),
            [&](int l) {
                return point_plane_signed_distance_gradient(
                           P[l], T0[l], T1[l], T2[l])
                    .eval();
            });
        check_lanes(
            "hessian", point_plane_signed_distance_hessian(p, t0, t1, t2),
            [&](int l) {
                return point_plane_signed_distance_hessian(
                           P[l], T0[l], T1[l], T2[l])
                    .eval();
            });
    }
}

#endif
