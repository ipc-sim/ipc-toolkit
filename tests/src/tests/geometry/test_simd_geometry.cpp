#include <catch2/catch_test_macros.hpp>

#include <tests/simd_utils.hpp>

#ifdef IPC_TOOLKIT_WITH_SIMD

#include <ipc/geometry/area.hpp>
#include <ipc/geometry/normal.hpp>

#include <cmath>
#include <string>

using namespace ipc;
using namespace ipc::tests;

TEST_CASE(
    "SIMD batch normals match the scalar ones lane-wise", "[normal][simd]")
{
    const Points<3> A = random_points<3>(1), B = random_points<3>(2),
                    C = random_points<3>(3), D = random_points<3>(4);

    const Eigen::Vector3<Batch> a = pack(A), b = pack(B), c = pack(C),
                                d = pack(D);

    SECTION("normalization")
    {
        check_lanes(
            "normalization_jacobian", normalization_jacobian(a),
            [&](int l) { return normalization_jacobian(A[l]).eval(); });
        // normalization_hessian returns one 3x3 slice per coordinate.
        const auto H = normalization_hessian(a);
        for (int k = 0; k < 3; ++k) {
            check_lanes(
                "normalization_hessian[" + std::to_string(k) + "]", H[k],
                [&](int l) { return normalization_hessian(A[l])[k]; });
        }
    }

    SECTION("point-line")
    {
        check_lanes(
            "unnormalized", point_line_unnormalized_normal(a, b, c),
            [&](int l) {
                return point_line_unnormalized_normal(A[l], B[l], C[l]).eval();
            });
        check_lanes("normal", point_line_normal(a, b, c), [&](int l) {
            return point_line_normal(A[l], B[l], C[l]).eval();
        });
        check_lanes(
            "unnormalized_jacobian",
            point_line_unnormalized_normal_jacobian(a, b, c), [&](int l) {
                return point_line_unnormalized_normal_jacobian(A[l], B[l], C[l])
                    .eval();
            });
        check_lanes(
            "normal_jacobian", point_line_normal_jacobian(a, b, c), [&](int l) {
                return point_line_normal_jacobian(A[l], B[l], C[l]).eval();
            });
        check_lanes(
            "unnormalized_hessian",
            point_line_unnormalized_normal_hessian(a, b, c), [&](int l) {
                return point_line_unnormalized_normal_hessian(A[l], B[l], C[l])
                    .eval();
            });
        check_lanes(
            "normal_hessian", point_line_normal_hessian(a, b, c), [&](int l) {
                return point_line_normal_hessian(A[l], B[l], C[l]).eval();
            });
    }

    SECTION("triangle")
    {
        check_lanes(
            "unnormalized", triangle_unnormalized_normal(a, b, c), [&](int l) {
                return triangle_unnormalized_normal(A[l], B[l], C[l]).eval();
            });
        check_lanes("normal", triangle_normal(a, b, c), [&](int l) {
            return triangle_normal(A[l], B[l], C[l]).eval();
        });
        check_lanes(
            "unnormalized_jacobian",
            triangle_unnormalized_normal_jacobian(a, b, c), [&](int l) {
                return triangle_unnormalized_normal_jacobian(A[l], B[l], C[l])
                    .eval();
            });
        check_lanes(
            "normal_jacobian", triangle_normal_jacobian(a, b, c), [&](int l) {
                return triangle_normal_jacobian(A[l], B[l], C[l]).eval();
            });
        check_lanes(
            "unnormalized_hessian",
            triangle_unnormalized_normal_hessian(a, b, c), [&](int l) {
                return triangle_unnormalized_normal_hessian(A[l], B[l], C[l])
                    .eval();
            });
        check_lanes(
            "normal_hessian", triangle_normal_hessian(a, b, c), [&](int l) {
                return triangle_normal_hessian(A[l], B[l], C[l]).eval();
            });
    }

    SECTION("line-line")
    {
        check_lanes(
            "unnormalized", line_line_unnormalized_normal(a, b, c, d),
            [&](int l) {
                return line_line_unnormalized_normal(A[l], B[l], C[l], D[l])
                    .eval();
            });
        check_lanes("normal", line_line_normal(a, b, c, d), [&](int l) {
            return line_line_normal(A[l], B[l], C[l], D[l]).eval();
        });
        check_lanes(
            "unnormalized_jacobian",
            line_line_unnormalized_normal_jacobian(a, b, c, d), [&](int l) {
                return line_line_unnormalized_normal_jacobian(
                           A[l], B[l], C[l], D[l])
                    .eval();
            });
        check_lanes(
            "normal_jacobian", line_line_normal_jacobian(a, b, c, d),
            [&](int l) {
                return line_line_normal_jacobian(A[l], B[l], C[l], D[l]).eval();
            });
        check_lanes(
            "unnormalized_hessian",
            line_line_unnormalized_normal_hessian(a, b, c, d), [&](int l) {
                return line_line_unnormalized_normal_hessian(
                           A[l], B[l], C[l], D[l])
                    .eval();
            });
        check_lanes(
            "normal_hessian", line_line_normal_hessian(a, b, c, d), [&](int l) {
                return line_line_normal_hessian(A[l], B[l], C[l], D[l]).eval();
            });
    }
}

TEST_CASE(
    "SIMD batch normalized normals leave a degenerate lane unscaled",
    "[normal][simd]")
{
    // ipc::normalized applies Eigen's own rule -- a zero-length vector is
    // returned unscaled -- per lane rather than behind one bool. Only the
    // values are compared: the normalized Jacobian/Hessian of a degenerate
    // configuration divides by zero on both paths.
    //
    // Lane 0 is degenerate and lane 1 is generic, and the roles swap on the
    // second pass, so both orderings of the blend are exercised even when a
    // batch holds only two lanes.
    for (int flip = 0; flip < 2; ++flip) {
        const auto is_degenerate = [&](int l) { return (l % 2) == flip; };

        // point-line: the point lies on the line, so d is parallel to e.
        // triangle: the three vertices are collinear.
        // line-line: the two lines are parallel.
        Points<3> P, E0, E1;
        Points<3> TA, TB, TC;
        Points<3> LA0, LA1, LB0, LB1;
        for (int l = 0; l < L; ++l) {
            // Exactly representable coordinates, so the degenerate lane's
            // unnormalized normal cancels to a true zero rather than to a
            // rounding residue that would take the non-degenerate branch.
            const Eigen::Vector3d e(1, 2, 3);
            E0[l] = Eigen::Vector3d(0.5 * l, -0.25, 0.75);
            E1[l] = E0[l] + e;
            P[l] = is_degenerate(l)
                ? Eigen::Vector3d(E0[l] + 0.5 * e)
                : Eigen::Vector3d(E0[l] + Eigen::Vector3d(0, 0, 1));

            TA[l] = Eigen::Vector3d(0, 0, 0);
            TB[l] = Eigen::Vector3d(1, 0, 0);
            TC[l] = is_degenerate(l) ? Eigen::Vector3d(2, 0, 0)
                                     : Eigen::Vector3d(0, 1, 0);

            LA0[l] = Eigen::Vector3d(0, 0, 0);
            LA1[l] = Eigen::Vector3d(1, 0, 0);
            LB0[l] = Eigen::Vector3d(0, 1, 0);
            LB1[l] = LB0[l]
                + (is_degenerate(l) ? Eigen::Vector3d(1, 0, 0)
                                    : Eigen::Vector3d(0, 0, 1));
        }

        check_lanes(
            "point_line_normal", point_line_normal(pack(P), pack(E0), pack(E1)),
            [&](int l) {
                return point_line_normal(P[l], E0[l], E1[l]).eval();
            });
        check_lanes(
            "triangle_normal", triangle_normal(pack(TA), pack(TB), pack(TC)),
            [&](int l) { return triangle_normal(TA[l], TB[l], TC[l]).eval(); });
        check_lanes(
            "line_line_normal",
            line_line_normal(pack(LA0), pack(LA1), pack(LB0), pack(LB1)),
            [&](int l) {
                return line_line_normal(LA0[l], LA1[l], LB0[l], LB1[l]).eval();
            });

        // Agreeing with the scalar path is not enough on its own: assert the
        // rule directly, so a degenerate lane is exactly zero rather than a
        // NaN, and a generic lane in the same batch is still a unit vector.
        const auto check_degenerate_rule = [&](const std::string& name,
                                               const Eigen::Vector3<Batch>& n) {
            for (int l = 0; l < L; ++l) {
                Eigen::Vector3d lane;
                for (int k = 0; k < 3; ++k) {
                    lane[k] = n[k].get(l);
                }
                CAPTURE(name, flip, l, lane.transpose());
                if (is_degenerate(l)) {
                    CHECK(lane == Eigen::Vector3d::Zero());
                } else {
                    CHECK(std::abs(lane.norm() - 1.0) <= DERIVATIVE_TOL);
                }
            }
        };

        check_degenerate_rule(
            "point_line_normal",
            point_line_normal(pack(P), pack(E0), pack(E1)));
        check_degenerate_rule(
            "triangle_normal", triangle_normal(pack(TA), pack(TB), pack(TC)));
        check_degenerate_rule(
            "line_line_normal",
            line_line_normal(pack(LA0), pack(LA1), pack(LB0), pack(LB1)));
    }
}

TEST_CASE(
    "SIMD batch edge length and triangle area match the scalar ones lane-wise",
    "[area][simd]")
{
    const Points<3> A = random_points<3>(5), B = random_points<3>(6),
                    C = random_points<3>(7);

    const Eigen::Vector3<Batch> a = pack(A), b = pack(B), c = pack(C);

    check_scalar_lanes("edge_length", edge_length(a, b), [&](int l) {
        return edge_length(A[l], B[l]);
    });
    check_scalar_lanes("triangle_area", triangle_area(a, b, c), [&](int l) {
        return triangle_area(A[l], B[l], C[l]);
    });

    check_lanes("edge_length_gradient", edge_length_gradient(a, b), [&](int l) {
        return edge_length_gradient(A[l], B[l]).eval();
    });
    check_lanes(
        "triangle_area_gradient", triangle_area_gradient(a, b, c),
        [&](int l) { return triangle_area_gradient(A[l], B[l], C[l]).eval(); });
}

#endif
