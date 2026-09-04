#include <catch2/catch_test_macros.hpp>

#include <tests/simd_utils.hpp>

#ifdef IPC_TOOLKIT_WITH_SIMD

#include <ipc/tangent/tangent_basis.hpp>

using namespace ipc;
using namespace ipc::tests;

TEST_CASE(
    "SIMD batch tangent bases match the scalar ones lane-wise",
    "[tangent_basis][simd]")
{
    const Points<3> A = random_points<3>(1), B = random_points<3>(2),
                    C = random_points<3>(3), D = random_points<3>(4);

    const Eigen::Vector3<Batch> a = pack(A), b = pack(B), c = pack(C),
                                d = pack(D);

    check_lanes("point_point", point_point_tangent_basis(a, b), [&](int l) {
        return point_point_tangent_basis(A[l], B[l]).eval();
    });
    check_lanes(
        "point_point_jacobian", point_point_tangent_basis_jacobian(a, b),
        [&](int l) {
            return point_point_tangent_basis_jacobian(A[l], B[l]).eval();
        });

    check_lanes("point_edge", point_edge_tangent_basis(a, b, c), [&](int l) {
        return point_edge_tangent_basis(A[l], B[l], C[l]).eval();
    });
    check_lanes(
        "point_edge_jacobian", point_edge_tangent_basis_jacobian(a, b, c),
        [&](int l) {
            return point_edge_tangent_basis_jacobian(A[l], B[l], C[l]).eval();
        });

    check_lanes("edge_edge", edge_edge_tangent_basis(a, b, c, d), [&](int l) {
        return edge_edge_tangent_basis(A[l], B[l], C[l], D[l]).eval();
    });
    check_lanes(
        "edge_edge_jacobian", edge_edge_tangent_basis_jacobian(a, b, c, d),
        [&](int l) {
            return edge_edge_tangent_basis_jacobian(A[l], B[l], C[l], D[l])
                .eval();
        });

    check_lanes(
        "point_triangle", point_triangle_tangent_basis(a, b, c, d), [&](int l) {
            return point_triangle_tangent_basis(A[l], B[l], C[l], D[l]).eval();
        });
    check_lanes(
        "point_triangle_jacobian",
        point_triangle_tangent_basis_jacobian(a, b, c, d), [&](int l) {
            return point_triangle_tangent_basis_jacobian(A[l], B[l], C[l], D[l])
                .eval();
        });
}

TEST_CASE(
    "SIMD batch point-point tangent basis blends the reference axis per lane",
    "[tangent_basis][simd]")
{
    // Lanes deliberately straddle the cross_x / cross_y branch: a pair along
    // x prefers one reference axis, a pair along y the other.
    Points<3> A, B;
    for (int l = 0; l < L; ++l) {
        A[l] = Eigen::Vector3d::Zero();
        B[l] =
            (l % 2 == 0) ? Eigen::Vector3d(1, 0, 0) : Eigen::Vector3d(0, 1, 0);
    }

    const Eigen::Vector3<Batch> a = pack(A), b = pack(B);

    check_lanes("point_point", point_point_tangent_basis(a, b), [&](int l) {
        return point_point_tangent_basis(A[l], B[l]).eval();
    });
    check_lanes(
        "point_point_jacobian", point_point_tangent_basis_jacobian(a, b),
        [&](int l) {
            return point_point_tangent_basis_jacobian(A[l], B[l]).eval();
        });
}

#endif
