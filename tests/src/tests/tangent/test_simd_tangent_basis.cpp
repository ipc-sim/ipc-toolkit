#include <catch2/catch_test_macros.hpp>

#include <ipc/utils/simd.hpp>

#ifdef IPC_TOOLKIT_WITH_SIMD

#include <ipc/tangent/tangent_basis.hpp>

#include <cmath>
#include <string>
#include <vector>

using namespace ipc;

namespace {

using Batch = SimdBatch<double>;
constexpr int L = int(Batch::size);

/// @brief The batch and scalar paths differ only in how the compiler contracts
/// multiply-adds. As with the distance derivatives, the tolerance is scaled by
/// the magnitude of the whole result rather than per-entry: a tangent-basis
/// Jacobian for a near-degenerate configuration has entries that are
/// catastrophic cancellations of much larger terms. A structurally wrong entry
/// would differ by order of the result's own magnitude.
constexpr double TOL = 1e-9;

/// @brief Pack the k-th coordinate of L problems into one batch.
Eigen::Vector3<Batch> pack(const std::vector<Eigen::Vector3d>& v)
{
    Eigen::Vector3<Batch> out;
    for (int k = 0; k < 3; ++k) {
        std::vector<double> tmp(L);
        for (int l = 0; l < L; ++l) {
            tmp[l] = v[l][k];
        }
        out[k] = Batch::load_unaligned(tmp.data());
    }
    return out;
}

std::vector<Eigen::Vector3d> random_points(const int seed)
{
    std::srand(seed);
    std::vector<Eigen::Vector3d> v(L);
    for (int l = 0; l < L; ++l) {
        v[l] = Eigen::Vector3d::Random();
    }
    return v;
}

/// @brief Compare a batch-valued matrix against the scalar result per lane.
template <typename BatchMatrix, typename ScalarOf>
void check_lanes(
    const std::string& name, const BatchMatrix& batched, ScalarOf&& scalar_of)
{
    for (int l = 0; l < L; ++l) {
        const auto want = scalar_of(l);
        REQUIRE(batched.rows() == want.rows());
        REQUIRE(batched.cols() == want.cols());
        const double scale = std::max(1.0, want.array().abs().maxCoeff());
        for (Eigen::Index i = 0; i < want.size(); ++i) {
            const double got = batched(i).get(l);
            CAPTURE(name, l, i, got, want(i), scale);
            CHECK(std::abs(got - want(i)) <= TOL * scale);
        }
    }
}

} // namespace

TEST_CASE(
    "SIMD batch tangent bases match the scalar ones lane-wise",
    "[tangent_basis][simd]")
{
    const std::vector<Eigen::Vector3d> A = random_points(1);
    const std::vector<Eigen::Vector3d> B = random_points(2);
    const std::vector<Eigen::Vector3d> C = random_points(3);
    const std::vector<Eigen::Vector3d> D = random_points(4);

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
    std::vector<Eigen::Vector3d> A(L, Eigen::Vector3d::Zero()), B(L);
    for (int l = 0; l < L; ++l) {
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
