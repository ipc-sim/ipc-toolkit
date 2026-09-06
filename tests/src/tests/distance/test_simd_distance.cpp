#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <tests/simd_utils.hpp>

#ifdef IPC_TOOLKIT_WITH_SIMD

#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/line_line.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/distance/point_line.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/distance/point_triangle.hpp>

using namespace ipc;
using namespace ipc::tests;

TEST_CASE(
    "SIMD batch distances match the scalar ones lane-wise", "[distance][simd]")
{
    const Points<3> A = random_points<3>(1), B = random_points<3>(2),
                    C = random_points<3>(3), D = random_points<3>(4);
    const Eigen::Vector3<Batch> a = pack(A), b = pack(B), c = pack(C),
                                d = pack(D);

    // A batch lane must agree with the scalar answer to within rounding; the
    // two differ only in how the compiler contracts multiply-adds.
    check_scalar_lanes(
        "point_point_distance", point_point_distance(a, b),
        [&](int l) { return point_point_distance(A[l], B[l]); });
    check_scalar_lanes(
        "point_line_distance", point_line_distance(a, b, c),
        [&](int l) { return point_line_distance(A[l], B[l], C[l]); });
    check_scalar_lanes(
        "line_line_distance", line_line_distance(a, b, c, d),
        [&](int l) { return line_line_distance(A[l], B[l], C[l], D[l]); });
    check_scalar_lanes(
        "point_edge_distance",
        point_edge_distance(a, b, c, PointEdgeDistanceType::P_E), [&](int l) {
            return point_edge_distance(
                A[l], B[l], C[l], PointEdgeDistanceType::P_E);
        });
    check_scalar_lanes(
        "edge_edge_distance",
        edge_edge_distance(a, b, c, d, EdgeEdgeDistanceType::EA_EB),
        [&](int l) {
            return edge_edge_distance(
                A[l], B[l], C[l], D[l], EdgeEdgeDistanceType::EA_EB);
        });
    check_scalar_lanes(
        "point_triangle_distance",
        point_triangle_distance(a, b, c, d, PointTriangleDistanceType::P_T),
        [&](int l) {
            return point_triangle_distance(
                A[l], B[l], C[l], D[l], PointTriangleDistanceType::P_T);
        });
}

TEST_CASE(
    "SIMD batch distance gradients and Hessians match lane-wise",
    "[distance][simd]")
{
    const Points<3> A = random_points<3>(11), B = random_points<3>(12),
                    C = random_points<3>(13), D = random_points<3>(14);
    const Eigen::Vector3<Batch> a = pack(A), b = pack(B), c = pack(C),
                                d = pack(D);

    check_lanes(
        "point_point_distance_gradient", point_point_distance_gradient(a, b),
        [&](int l) {
            return point_point_distance_gradient(A[l], B[l]).eval();
        });
    check_lanes(
        "point_point_distance_hessian", point_point_distance_hessian(a, b),
        [&](int l) { return point_point_distance_hessian(A[l], B[l]).eval(); });

    check_lanes(
        "point_line_distance_gradient", point_line_distance_gradient(a, b, c),
        [&](int l) {
            return point_line_distance_gradient(A[l], B[l], C[l]).eval();
        });
    check_lanes(
        "point_line_distance_hessian", point_line_distance_hessian(a, b, c),
        [&](int l) {
            return point_line_distance_hessian(A[l], B[l], C[l]).eval();
        });

    check_lanes(
        "line_line_distance_gradient", line_line_distance_gradient(a, b, c, d),
        [&](int l) {
            return line_line_distance_gradient(A[l], B[l], C[l], D[l]).eval();
        });
    check_lanes(
        "line_line_distance_hessian", line_line_distance_hessian(a, b, c, d),
        [&](int l) {
            return line_line_distance_hessian(A[l], B[l], C[l], D[l]).eval();
        });

    check_lanes(
        "point_edge_distance_gradient",
        point_edge_distance_gradient(a, b, c, PointEdgeDistanceType::P_E),
        [&](int l) {
            return point_edge_distance_gradient(
                       A[l], B[l], C[l], PointEdgeDistanceType::P_E)
                .eval();
        });
    check_lanes(
        "point_edge_distance_hessian",
        point_edge_distance_hessian(a, b, c, PointEdgeDistanceType::P_E),
        [&](int l) {
            return point_edge_distance_hessian(
                       A[l], B[l], C[l], PointEdgeDistanceType::P_E)
                .eval();
        });

    check_lanes(
        "edge_edge_distance_gradient",
        edge_edge_distance_gradient(a, b, c, d, EdgeEdgeDistanceType::EA_EB),
        [&](int l) {
            return edge_edge_distance_gradient(
                       A[l], B[l], C[l], D[l], EdgeEdgeDistanceType::EA_EB)
                .eval();
        });
    check_lanes(
        "edge_edge_distance_hessian",
        edge_edge_distance_hessian(a, b, c, d, EdgeEdgeDistanceType::EA_EB),
        [&](int l) {
            return edge_edge_distance_hessian(
                       A[l], B[l], C[l], D[l], EdgeEdgeDistanceType::EA_EB)
                .eval();
        });

    check_lanes(
        "point_triangle_distance_gradient",
        point_triangle_distance_gradient(
            a, b, c, d, PointTriangleDistanceType::P_T),
        [&](int l) {
            return point_triangle_distance_gradient(
                       A[l], B[l], C[l], D[l], PointTriangleDistanceType::P_T)
                .eval();
        });
    check_lanes(
        "point_triangle_distance_hessian",
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
    const Points<3> A = random_points<3>(5), B = random_points<3>(6),
                    C = random_points<3>(7), D = random_points<3>(8);
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
