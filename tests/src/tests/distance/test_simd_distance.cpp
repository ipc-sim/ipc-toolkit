#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <ipc/utils/simd.hpp>

#ifdef IPC_TOOLKIT_WITH_SIMD

#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/line_line.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/distance/point_line.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/distance/point_triangle.hpp>

#include <cmath>
#include <vector>

using namespace ipc;

namespace {

using Batch = SimdBatch<double>;
constexpr int L = int(Batch::size);

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

/// @brief Agreement to within rounding; the batch and scalar paths differ only
/// in how the compiler contracts multiply-adds.
bool close(const double got, const double want)
{
    return std::abs(got - want) <= 1e-14 * std::max(1.0, std::abs(want));
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

} // namespace

TEST_CASE(
    "SIMD batch distances match the scalar ones lane-wise", "[distance][simd]")
{
    const std::vector<Eigen::Vector3d> A = random_points(1),
                                       B = random_points(2),
                                       C = random_points(3),
                                       D = random_points(4);
    const Eigen::Vector3<Batch> a = pack(A), b = pack(B), c = pack(C),
                                d = pack(D);

    const Batch pp = point_point_distance(a, b);
    const Batch pl = point_line_distance(a, b, c);
    const Batch ll = line_line_distance(a, b, c, d);
    const Batch pe = point_edge_distance(a, b, c, PointEdgeDistanceType::P_E);
    const Batch ee =
        edge_edge_distance(a, b, c, d, EdgeEdgeDistanceType::EA_EB);
    const Batch pt =
        point_triangle_distance(a, b, c, d, PointTriangleDistanceType::P_T);

    for (int l = 0; l < L; ++l) {
        CAPTURE(l);
        // A batch lane must agree with the scalar answer to within rounding;
        // the two differ only in how the compiler contracts multiply-adds.
        CHECK(close(pp.get(l), point_point_distance(A[l], B[l])));
        CHECK(close(pl.get(l), point_line_distance(A[l], B[l], C[l])));
        CHECK(close(ll.get(l), line_line_distance(A[l], B[l], C[l], D[l])));
        CHECK(close(
            pe.get(l),
            point_edge_distance(A[l], B[l], C[l], PointEdgeDistanceType::P_E)));
        CHECK(close(
            ee.get(l),
            edge_edge_distance(
                A[l], B[l], C[l], D[l], EdgeEdgeDistanceType::EA_EB)));
        CHECK(close(
            pt.get(l),
            point_triangle_distance(
                A[l], B[l], C[l], D[l], PointTriangleDistanceType::P_T)));
    }
}

TEST_CASE(
    "SIMD batch distance gradients and Hessians match lane-wise",
    "[distance][simd]")
{
    const std::vector<Eigen::Vector3d> A = random_points(11),
                                       B = random_points(12),
                                       C = random_points(13),
                                       D = random_points(14);
    const Eigen::Vector3<Batch> a = pack(A), b = pack(B), c = pack(C),
                                d = pack(D);

    // Compare a batched derivative against the scalar one, lane by lane.
    //
    // NOTE: the tolerance is scaled by the magnitude of the whole result, not
    // by each entry. These derivatives are ill-conditioned for near-parallel
    // edges -- a random configuration routinely produces entries of order 1e6
    // whose small differences are catastrophic cancellations of large terms.
    // The batch and scalar paths contract multiply-adds differently, so they
    // agree to roughly 1e-11 relative to the norm rather than to an ulp. A
    // structurally wrong entry would differ by order of the norm itself.
    const auto check_all = [&](const auto& batched, const auto& scalar_of) {
        for (int l = 0; l < L; ++l) {
            CAPTURE(l);
            const auto want = scalar_of(l);
            REQUIRE(batched.size() == want.size());
            const double scale = std::max(1.0, want.array().abs().maxCoeff());
            for (Eigen::Index i = 0; i < want.size(); ++i) {
                CAPTURE(i);
                const double got = batched(i).get(l), expect = want(i);
                CAPTURE(got, expect, std::abs(got - expect), scale);
                CHECK(std::abs(got - expect) <= 1e-9 * scale);
            }
        }
    };

    check_all(point_point_distance_gradient(a, b), [&](int l) {
        return point_point_distance_gradient(A[l], B[l]).eval();
    });
    check_all(point_point_distance_hessian(a, b), [&](int l) {
        return point_point_distance_hessian(A[l], B[l]).eval();
    });

    check_all(point_line_distance_gradient(a, b, c), [&](int l) {
        return point_line_distance_gradient(A[l], B[l], C[l]).eval();
    });
    check_all(point_line_distance_hessian(a, b, c), [&](int l) {
        return point_line_distance_hessian(A[l], B[l], C[l]).eval();
    });

    check_all(line_line_distance_gradient(a, b, c, d), [&](int l) {
        return line_line_distance_gradient(A[l], B[l], C[l], D[l]).eval();
    });
    check_all(line_line_distance_hessian(a, b, c, d), [&](int l) {
        return line_line_distance_hessian(A[l], B[l], C[l], D[l]).eval();
    });

    check_all(
        point_edge_distance_gradient(a, b, c, PointEdgeDistanceType::P_E),
        [&](int l) {
            return point_edge_distance_gradient(
                       A[l], B[l], C[l], PointEdgeDistanceType::P_E)
                .eval();
        });
    check_all(
        point_edge_distance_hessian(a, b, c, PointEdgeDistanceType::P_E),
        [&](int l) {
            return point_edge_distance_hessian(
                       A[l], B[l], C[l], PointEdgeDistanceType::P_E)
                .eval();
        });

    check_all(
        edge_edge_distance_gradient(a, b, c, d, EdgeEdgeDistanceType::EA_EB),
        [&](int l) {
            return edge_edge_distance_gradient(
                       A[l], B[l], C[l], D[l], EdgeEdgeDistanceType::EA_EB)
                .eval();
        });
    check_all(
        edge_edge_distance_hessian(a, b, c, d, EdgeEdgeDistanceType::EA_EB),
        [&](int l) {
            return edge_edge_distance_hessian(
                       A[l], B[l], C[l], D[l], EdgeEdgeDistanceType::EA_EB)
                .eval();
        });

    check_all(
        point_triangle_distance_gradient(
            a, b, c, d, PointTriangleDistanceType::P_T),
        [&](int l) {
            return point_triangle_distance_gradient(
                       A[l], B[l], C[l], D[l], PointTriangleDistanceType::P_T)
                .eval();
        });
    check_all(
        point_triangle_distance_hessian(
            a, b, c, d, PointTriangleDistanceType::P_T),
        [&](int l) {
            return point_triangle_distance_hessian(
                       A[l], B[l], C[l], D[l], PointTriangleDistanceType::P_T)
                .eval();
        });
}

TEST_CASE("SIMD batch distances reject AUTO", "[distance][simd]")
{
    // The distance type is a per-lane property, but the predicates return a
    // single enum, so AUTO cannot be resolved for a batch. It must throw
    // rather than silently apply one lane's classification to all of them.
    const std::vector<Eigen::Vector3d> A = random_points(5),
                                       B = random_points(6),
                                       C = random_points(7),
                                       D = random_points(8);
    const Eigen::Vector3<Batch> a = pack(A), b = pack(B), c = pack(C),
                                d = pack(D);

    // Not just any exception: the specific diagnostic naming the function,
    // rather than the generic "invalid distance type" from the switch default.
    const auto explicit_required = Catch::Matchers::MessageMatches(
        Catch::Matchers::ContainsSubstring("explicit distance type"));

    CHECK_THROWS_MATCHES(
        point_edge_distance(a, b, c, PointEdgeDistanceType::AUTO),
        std::invalid_argument, explicit_required);
    CHECK_THROWS_MATCHES(
        edge_edge_distance(a, b, c, d, EdgeEdgeDistanceType::AUTO),
        std::invalid_argument, explicit_required);
    CHECK_THROWS_MATCHES(
        point_triangle_distance(a, b, c, d, PointTriangleDistanceType::AUTO),
        std::invalid_argument, explicit_required);
    CHECK_THROWS_MATCHES(
        edge_edge_distance_gradient(a, b, c, d, EdgeEdgeDistanceType::AUTO),
        std::invalid_argument, explicit_required);
    CHECK_THROWS_MATCHES(
        edge_edge_distance_hessian(a, b, c, d, EdgeEdgeDistanceType::AUTO),
        std::invalid_argument, explicit_required);
    CHECK_THROWS_MATCHES(
        point_triangle_distance_gradient(
            a, b, c, d, PointTriangleDistanceType::AUTO),
        std::invalid_argument, explicit_required);
    CHECK_THROWS_MATCHES(
        point_triangle_distance_hessian(
            a, b, c, d, PointTriangleDistanceType::AUTO),
        std::invalid_argument, explicit_required);
}

#endif
