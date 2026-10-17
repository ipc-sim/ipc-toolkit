#include <tests/config.hpp>
#include <tests/utils.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <ipc/esp/arbitrary_point_esp.hpp>

#include <finitediff.hpp>
#include <igl/edges.h>

#include <algorithm>
#include <cmath>

using namespace ipc;

namespace {

struct Fixture {
    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    CollisionMesh mesh;
    double dhat;

    Fixture()
    {
        REQUIRE(tests::load_mesh("cube.ply", V, E, F));
        mesh = CollisionMesh(V, E, F);
        const double bbox_diag =
            (V.colwise().maxCoeff() - V.colwise().minCoeff()).norm();
        dhat = 0.2 * bbox_diag;
    }

    // A point just off the surface, offset from face 0's centroid along its
    // (unnormalized-triangle) outward normal by a fraction of dhat.
    Eigen::RowVector3d near_surface_point(double frac = 0.3) const
    {
        const Eigen::RowVector3d v0 = V.row(F(0, 0));
        const Eigen::RowVector3d v1 = V.row(F(0, 1));
        const Eigen::RowVector3d v2 = V.row(F(0, 2));
        const Eigen::RowVector3d centroid = (v0 + v1 + v2) / 3.0;
        const Eigen::RowVector3d normal = (v1 - v0).cross(v2 - v0).normalized();
        return centroid + frac * dhat * normal;
    }
};

// The 2D analogue: the unit square as a closed 4-edge polyline, so a query
// point can be placed against an edge interior or against a corner (where two
// edges both reduce to the same vertex and must cancel against the direct
// vertex term).
struct Fixture2D {
    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    CollisionMesh mesh;
    double dhat = 0.2;

    Fixture2D()
    {
        V.resize(4, 2);
        V << 0, 0, 1, 0, 1, 1, 0, 1;
        E.resize(4, 2);
        E << 0, 1, 1, 2, 2, 3, 3, 0;
        mesh = CollisionMesh(V, E, F);
    }

    // Outside the bottom edge's midpoint (closest feature is an edge
    // interior), at frac * dhat.
    Eigen::RowVector2d near_edge_point(double frac) const
    {
        return Eigen::RowVector2d(0.5, -frac * dhat);
    }

    // Outside the corner at vertex 0, along the diagonal (closest feature is
    // a corner shared by two edges).
    Eigen::RowVector2d near_corner_point(double frac) const
    {
        return Eigen::RowVector2d(0, 0)
            - frac * dhat * Eigen::RowVector2d(1, 1).normalized();
    }
};

} // namespace

TEST_CASE(
    "Arbitrary Point ESP: zero beyond dhat",
    "[esp_potential],[arbitrary_point_esp]")
{
    Fixture fx;
    ESPParameters params(fx.dhat);
    ArbitraryPointESP<3> potential(fx.mesh, params);
    potential.update(fx.V);

    const Eigen::RowVector3d far_point =
        fx.V.colwise().maxCoeff() + Eigen::RowVector3d::Constant(10 * fx.dhat);

    REQUIRE(potential(fx.V, far_point) == 0.0);
    REQUIRE(potential.gradient(fx.V, far_point).isZero());
    REQUIRE(potential.hessian(fx.V, far_point).isZero());
}

TEST_CASE(
    "Arbitrary Point ESP: FD gradient/hessian at an off-mesh point",
    "[esp_potential],[arbitrary_point_esp]")
{
    Fixture fx;
    ESPParameters params(fx.dhat);
    ArbitraryPointESP<3> potential(fx.mesh, params);
    potential.update(fx.V);

    const Eigen::RowVector3d q = fx.near_surface_point();

    // Sanity: the point should actually be within range of the surface.
    REQUIRE(potential(fx.V, q) != 0.0);

    SECTION("gradient")
    {
        const Eigen::Vector3d g = potential.gradient(fx.V, q);

        Eigen::VectorXd fg;
        fd::finite_gradient(
            Eigen::VectorXd(q.transpose()),
            [&](const Eigen::VectorXd& y) {
                return potential(fx.V, Eigen::RowVector3d(y.transpose()));
            },
            fg, fd::AccuracyOrder::SECOND, 1e-8);

        REQUIRE((fg - g).norm() < std::max(1e-8, fg.norm()) * 1e-5);
    }

    SECTION("hessian")
    {
        const Eigen::Matrix3d h = potential.hessian(fx.V, q);

        Eigen::MatrixXd fh;
        fd::finite_jacobian(
            Eigen::VectorXd(q.transpose()),
            [&](const Eigen::VectorXd& y) {
                return potential.gradient(
                    fx.V, Eigen::RowVector3d(y.transpose()));
            },
            fh, fd::AccuracyOrder::SECOND, 1e-8);

        REQUIRE((fh - h).norm() < std::max(1e-8, fh.norm()) * 1e-5);
    }
}

TEST_CASE(
    "Arbitrary Point ESP: evaluate() matches operator()/gradient()/hessian()",
    "[esp_potential],[arbitrary_point_esp]")
{
    Fixture fx;
    ESPParameters params(fx.dhat);
    ArbitraryPointESP<3> potential(fx.mesh, params);
    potential.update(fx.V);

    // Sweep from just outside dhat down through the surface to well inside
    // the mesh, so both the "no collisions" and "several collisions merged
    // via symbolic cancellation" code paths are exercised.
    for (const double frac : { -0.3, 0.05, 0.3, 0.7, 0.95 }) {
        CAPTURE(frac);
        const Eigen::RowVector3d q = fx.near_surface_point(frac);

        const double value = potential(fx.V, q);
        const Eigen::Vector3d grad = potential.gradient(fx.V, q);
        const Eigen::Matrix3d hess = potential.hessian(fx.V, q);

        const auto [value2, grad2, hess2] = potential.evaluate(fx.V, q);

        // evaluate() runs the exact same per-collision computation and
        // accumulation order as the three separate calls, just fused into
        // one pass -- expect bit-exact agreement, not just Approx.
        REQUIRE(value2 == value);
        REQUIRE(grad2 == grad);
        REQUIRE(hess2 == hess);
    }

    // Beyond dhat: all three outputs zero, consistent with the separate
    // accessors (see the "zero beyond dhat" test above).
    const Eigen::RowVector3d far_point =
        fx.V.colwise().maxCoeff() + Eigen::RowVector3d::Constant(10 * fx.dhat);
    const auto [value, grad, hess] = potential.evaluate(fx.V, far_point);
    REQUIRE(value == 0.0);
    REQUIRE(grad.isZero());
    REQUIRE(hess.isZero());
}

TEST_CASE(
    "Arbitrary Point ESP 2D: zero beyond dhat",
    "[esp_potential],[arbitrary_point_esp]")
{
    Fixture2D fx;
    ESPParameters params(fx.dhat);
    ArbitraryPointESP<2> potential(fx.mesh, params);
    potential.update(fx.V);

    const Eigen::RowVector2d far_point(0.5, -10 * fx.dhat);

    REQUIRE(potential(fx.V, far_point) == 0.0);
    REQUIRE(potential.gradient(fx.V, far_point).isZero());
    REQUIRE(potential.hessian(fx.V, far_point).isZero());
}

TEST_CASE(
    "Arbitrary Point ESP 2D: FD gradient/hessian at an off-mesh point",
    "[esp_potential],[arbitrary_point_esp]")
{
    Fixture2D fx;
    ESPParameters params(fx.dhat);
    ArbitraryPointESP<2> potential(fx.mesh, params);
    potential.update(fx.V);

    // Both an edge-interior closest feature and a corner, where the two
    // incident edges reduce to the same vertex and their +1s must merge with
    // the direct vertex term's -1.
    const Eigen::RowVector2d q = GENERATE_COPY(
        Eigen::RowVector2d(fx.near_edge_point(0.3)),
        Eigen::RowVector2d(fx.near_corner_point(0.3)));
    CAPTURE(q);

    REQUIRE(potential(fx.V, q) != 0.0);

    SECTION("gradient")
    {
        const Eigen::Vector2d g = potential.gradient(fx.V, q);

        Eigen::VectorXd fg;
        fd::finite_gradient(
            Eigen::VectorXd(q.transpose()),
            [&](const Eigen::VectorXd& y) {
                return potential(fx.V, Eigen::RowVector2d(y.transpose()));
            },
            fg, fd::AccuracyOrder::SECOND, 1e-8);

        REQUIRE((fg - g).norm() < std::max(1e-8, fg.norm()) * 1e-5);
    }

    SECTION("hessian")
    {
        const Eigen::Matrix2d h = potential.hessian(fx.V, q);

        Eigen::MatrixXd fh;
        fd::finite_jacobian(
            Eigen::VectorXd(q.transpose()),
            [&](const Eigen::VectorXd& y) {
                return potential.gradient(
                    fx.V, Eigen::RowVector2d(y.transpose()));
            },
            fh, fd::AccuracyOrder::SECOND, 1e-8);

        REQUIRE((fh - h).norm() < std::max(1e-8, fh.norm()) * 1e-5);
    }
}

TEST_CASE(
    "Arbitrary Point ESP 2D: corner value is a single vertex-vertex term",
    "[esp_potential],[arbitrary_point_esp]")
{
    // Outside the convex corner at V.row(0), both incident edges reduce to
    // that corner vertex (+1 each) and the direct vertex term contributes -1,
    // so the symbolic cancellation must leave exactly one vertex-vertex
    // barrier at the corner distance. This is the assertion that fails if the
    // 2D sign convention or the endpoint reduction is wrong: getting it
    // backwards leaves 3 terms or 0, both of which still pass a finite
    // difference check.
    Fixture2D fx;
    ESPParameters params(fx.dhat);
    ArbitraryPointESP<2> potential(fx.mesh, params);
    potential.update(fx.V);

    for (const double frac : { 0.2, 0.5, 0.9 }) {
        CAPTURE(frac);
        const Eigen::RowVector2d q = fx.near_corner_point(frac);
        const double d = (q - fx.V.row(0)).norm();
        REQUIRE(
            potential(fx.V, q) == Catch::Approx((*params.barrier)(d, fx.dhat)));
    }
}

TEST_CASE(
    "Arbitrary Point ESP: inclusion-exclusion weights on non-closed inputs",
    "[esp_potential],[arbitrary_point_esp]")
{
    // Inputs that are not closed manifolds: an open sheet and an open
    // polyline (boundary edges, boundary vertices and polyline ends weigh
    // zero, supplemental S4), a 3D polyline in no face (its edges weigh one,
    // its interior vertices minus one), isolated vertices (one) and a 2D
    // junction of three edges (minus two). On each input the query points lie
    // where the potential is the barrier of the distance. With the
    // closed-mesh signs instead, the potential is exactly zero beyond a
    // boundary edge, and edges in no face and 2D isolated vertices enter as
    // negative barriers.
    const double dhat = 0.5;
    ESPParameters params(dhat);

    SECTION("3D sheet")
    {
        // The unit square at z = 0 as a 3x3 grid of vertices, 8 triangles.
        Eigen::MatrixXd V(9, 3);
        for (int i = 0; i < 9; i++) {
            V.row(i) << (i % 3) / 2.0, (i / 3) / 2.0, 0;
        }
        Eigen::MatrixXi F(8, 3);
        F << 0, 1, 4, 0, 4, 3, 1, 2, 5, 1, 5, 4, 3, 4, 7, 3, 7, 6, 4, 5, 8, 4,
            8, 7;
        Eigen::MatrixXi E;
        igl::edges(F, E);
        const CollisionMesh mesh(V, E, F);
        ArbitraryPointESP<3> potential(mesh, params);
        potential.update(V);

        const Eigen::RowVector3d q = GENERATE(
            Eigen::RowVector3d(1.2, 0.5, 0.0),  // beyond a boundary edge
            Eigen::RowVector3d(1.1, 0.5, 0.2),  // beyond it, above the plane
            Eigen::RowVector3d(1.2, 1.1, 0.1),  // beyond a boundary corner
            Eigen::RowVector3d(0.9, 0.5, 0.1),  // inside, near the boundary
            Eigen::RowVector3d(1.0, 0.5, 0.2),  // above a boundary vertex
            Eigen::RowVector3d(0.5, 0.5, 0.2)); // above the interior vertex
        CAPTURE(q);

        const double dx = std::max({ 0.0, -q(0), q(0) - 1 });
        const double dy = std::max({ 0.0, -q(1), q(1) - 1 });
        const double d = std::sqrt(dx * dx + dy * dy + q(2) * q(2));
        CHECK(potential(V, q) == Catch::Approx((*params.barrier)(d, dhat)));

        Eigen::VectorXd fg;
        fd::finite_gradient(
            Eigen::VectorXd(q.transpose()),
            [&](const Eigen::VectorXd& y) {
                return potential(V, Eigen::RowVector3d(y.transpose()));
            },
            fg, fd::AccuracyOrder::SECOND, 1e-8);
        const Eigen::Vector3d g = potential.gradient(V, q);
        CHECK((fg - g).norm() < std::max(1e-8, fg.norm()) * 1e-5);
    }

    SECTION("2D polyline")
    {
        // The segment [0, 1] x {0} as three edges.
        Eigen::MatrixXd V(4, 2);
        V << 0, 0, 1 / 3.0, 0, 2 / 3.0, 0, 1, 0;
        Eigen::MatrixXi E(3, 2);
        E << 0, 1, 1, 2, 2, 3;
        const CollisionMesh mesh(V, E, Eigen::MatrixXi());
        ArbitraryPointESP<2> potential(mesh, params);
        potential.update(V);

        const Eigen::RowVector2d q = GENERATE(
            Eigen::RowVector2d(-0.2, 0.0),    // beyond an end
            Eigen::RowVector2d(-0.1, 0.1),    // beyond it, off the line
            Eigen::RowVector2d(0.1, 0.2),     // beside the line, near an end
            Eigen::RowVector2d(1 / 3.0, 0.1), // above an interior vertex
            Eigen::RowVector2d(0.5, -0.15));  // beside an edge interior
        CAPTURE(q);

        const double dx = std::max({ 0.0, -q(0), q(0) - 1 });
        const double d = std::sqrt(dx * dx + q(1) * q(1));
        CHECK(potential(V, q) == Catch::Approx((*params.barrier)(d, dhat)));
    }

    SECTION("3D polyline in no face")
    {
        // The segment [0, 1] x {0} x {0} as three edges, no faces.
        Eigen::MatrixXd V(4, 3);
        V << 0, 0, 0, 1 / 3.0, 0, 0, 2 / 3.0, 0, 0, 1, 0, 0;
        Eigen::MatrixXi E(3, 2);
        E << 0, 1, 1, 2, 2, 3;
        const CollisionMesh mesh(V, E, Eigen::MatrixXi());
        ArbitraryPointESP<3> potential(mesh, params);
        potential.update(V);

        const Eigen::RowVector3d q = GENERATE(
            Eigen::RowVector3d(-0.2, 0.0, 0.0),    // beyond an end
            Eigen::RowVector3d(0.1, 0.1, 0.1),     // beside it, near an end
            Eigen::RowVector3d(1 / 3.0, 0.0, 0.2), // beside an interior vertex
            Eigen::RowVector3d(0.5, -0.1, 0.1));   // beside an edge interior
        CAPTURE(q);

        const double dx = std::max({ 0.0, -q(0), q(0) - 1 });
        const double d = std::sqrt(dx * dx + q(1) * q(1) + q(2) * q(2));
        CHECK(potential(V, q) == Catch::Approx((*params.barrier)(d, dhat)));
    }

    SECTION("isolated vertices")
    {
        const Eigen::MatrixXd V3 = Eigen::MatrixXd::Zero(1, 3);
        const CollisionMesh mesh3(V3, Eigen::MatrixXi(), Eigen::MatrixXi());
        ArbitraryPointESP<3> potential3(mesh3, params);
        potential3.update(V3);
        const Eigen::RowVector3d q3(0.1, -0.2, 0.15);
        CHECK(
            potential3(V3, q3)
            == Catch::Approx((*params.barrier)(q3.norm(), dhat)));

        const Eigen::MatrixXd V2 = Eigen::MatrixXd::Zero(1, 2);
        const CollisionMesh mesh2(V2, Eigen::MatrixXi(), Eigen::MatrixXi());
        ArbitraryPointESP<2> potential2(mesh2, params);
        potential2.update(V2);
        const Eigen::RowVector2d q2(0.1, -0.2);
        CHECK(
            potential2(V2, q2)
            == Catch::Approx((*params.barrier)(q2.norm(), dhat)));
    }

    SECTION("2D junction")
    {
        // A bar (-1, 0)-(0, 0)-(1, 0) and a stem (0, 0)-(0, 1). Below the bar
        // the ball around a query point meets the input in one connected
        // piece, so the potential is the barrier of the distance.
        Eigen::MatrixXd V(4, 2);
        V << -1, 0, 0, 0, 1, 0, 0, 1;
        Eigen::MatrixXi E(3, 2);
        E << 0, 1, 1, 2, 1, 3;
        const CollisionMesh mesh(V, E, Eigen::MatrixXi());
        ArbitraryPointESP<2> potential(mesh, params);
        potential.update(V);

        const Eigen::RowVector2d q = GENERATE(
            Eigen::RowVector2d(0.0, -0.2),    // below the junction
            Eigen::RowVector2d(0.1, -0.1),    // below, beside the junction
            Eigen::RowVector2d(-0.3, -0.25)); // below the bar
        CAPTURE(q);

        CHECK(potential(V, q) == Catch::Approx((*params.barrier)(-q(1), dhat)));
    }
}

TEST_CASE(
    "Arbitrary Point ESP: mesh dimension mismatch",
    "[esp_potential][arbitrary_point_esp]")
{
    const Fixture2D fx2;
    CHECK_THROWS_AS(
        ArbitraryPointESP<3>(fx2.mesh, ESPParameters(fx2.dhat)),
        std::runtime_error);
    const Fixture fx3;
    CHECK_THROWS_AS(
        ArbitraryPointESP<2>(fx3.mesh, ESPParameters(fx3.dhat)),
        std::runtime_error);
}
