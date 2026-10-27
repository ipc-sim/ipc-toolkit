#include <tests/dof_layout.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/catch_approx.hpp>
#include <finitediff.hpp>
#include <igl/edges.h>
#include <ipc/esp/esp_potential.hpp>
#include <ipc/collisions/tangential/tangential_collisions.hpp>
#include <Eigen/Eigenvalues>

using namespace ipc;

namespace {
struct Scene {
    Eigen::MatrixXd rest, X;
    Eigen::MatrixXi E, F;
    CollisionMesh mesh;
};
Scene edge_over_box()
{
    Scene s;
    s.rest.resize(10, 3);
    s.rest << .035, -.016, .0006, .065, .017, .0008, -.05, -.04, -.1, .05, -.04,
        -.1, .05, .04, -.1, -.05, .04, -.1, -.05, -.04, 0, .05, -.04, 0, .05,
        .04, 0, -.05, .04, 0;
    s.F.resize(12, 3);
    s.F << 0, 3, 2, 0, 2, 1, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 3, 7, 6, 3, 6,
        2, 0, 4, 7, 0, 7, 3, 1, 2, 6, 1, 6, 5;
    s.F.array() += 2;
    Eigen::MatrixXi box_edges;
    igl::edges(s.F, box_edges);
    s.E.resize(box_edges.rows() + 1, 2);
    s.E.row(0) << 0, 1;
    s.E.bottomRows(box_edges.rows()) = box_edges;
    s.mesh = CollisionMesh(s.rest, s.E, s.F);
    Eigen::VectorXi patches(10);
    patches << 0, 0, 1, 1, 1, 1, 1, 1, 1, 1;
    s.mesh.can_collide = make_vertex_patches_filter(patches);
    s.mesh.init_adjacencies();
    s.X = s.rest;
    s.X(0, 0) = .0499615;
    return s;
}

void check_derivatives(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& X,
    const ESPParameters& params)
{
    const ESPPotential potential(params);
    ESPCollisions c;
    c.build(mesh, X, params);
    const auto g = potential.gradient(c, mesh, X);
    const Eigen::MatrixXd H = potential.hessian(c, mesh, X);
    Eigen::VectorXd fg;
    fd::finite_gradient(
        tests::flatten(X),
        [&](const Eigen::VectorXd& x) {
            const Eigen::MatrixXd V = tests::unflatten(x, 3);
            ESPCollisions contacts;
            contacts.build(mesh, V, params);
            return potential(contacts, mesh, V);
        },
        fg, fd::AccuracyOrder::FOURTH, 1e-8);
    CHECK((fg - g).norm() < 2e-5 * std::max(1.0, g.norm()));
    Eigen::MatrixXd fH;
    fd::finite_jacobian(
        tests::flatten(X),
        [&](const Eigen::VectorXd& x) {
            const Eigen::MatrixXd V = tests::unflatten(x, 3);
            ESPCollisions contacts;
            contacts.build(mesh, V, params);
            return potential.gradient(contacts, mesh, V);
        },
        fH, fd::AccuracyOrder::FOURTH, 1e-8);
    CHECK((fH - H).norm() < 2e-5 * std::max(1.0, H.norm()));
    CHECK((H - H.transpose()).norm() < 1e-10 * std::max(1.0, H.norm()));
    const Eigen::MatrixXd psd =
        potential.hessian(c, mesh, X, PSDProjectionMethod::CLAMP);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(psd);
    REQUIRE(eig.info() == Eigen::Success);
    CHECK(eig.eigenvalues().minCoeff() >= -1e-9 * std::max(1.0, psd.norm()));

    // Translation invariance catches missing endpoint/reference derivatives.
    const Eigen::MatrixXd forces = tests::unflatten(g, 3);
    CHECK(forces.colwise().sum().norm() < 1e-10 * std::max(1.0, g.norm()));
}
} // namespace

TEST_CASE(
    "ESP local transfer matches the edge-over-box prototype",
    "[esp_potential][esp_local_transfer]")
{
    auto s = edge_over_box();
    ESPParameters old_params(.03, .05);
    ESPParameters params = old_params;
    params.local_quadrature_transfer = true;
    const ESPPotential old_p(old_params), new_p(params);
    ESPCollisions c;
    c.build(s.mesh, s.X, params);
    const int ax = tests::vertex_dofs<3>(0, s.X.rows())[0];
    CHECK(
        old_p(c, s.mesh, s.X)
        == Catch::Approx(.10589380670845851).epsilon(1e-9));
    CHECK(
        new_p(c, s.mesh, s.X)
        == Catch::Approx(.0893123239682775).epsilon(1e-9));
    CHECK(
        -old_p.gradient(c, s.mesh, s.X)[ax]
        == Catch::Approx(852.69739536).epsilon(1e-8));
    CHECK(
        -new_p.gradient(c, s.mesh, s.X)[ax]
        == Catch::Approx(77.01422187).epsilon(1e-7));
}

TEST_CASE(
    "ESP local transfer full gradient Hessian and PSD",
    "[esp_potential][esp_local_transfer]")
{
    auto s = edge_over_box();
    ESPParameters params(.03, .05);
    params.local_quadrature_transfer = true;
    params.area_weights = GENERATE(true, false);
    const int barrier = GENERATE(0, 1, 2);
    if (barrier)
        params.barrier = std::make_shared<InversePowerBarrier<>>(barrier);
    // The vertex rule, corners plus a centroid and edge midpoints, or
    // per-face lattice rules: the corners donate, the other points are
    // ordinary samples.
    const int rule_kind = GENERATE(0, 1, 2);
    const bool with_non_corner_points = rule_kind == 1;
    if (rule_kind == 2) {
        params.set_auto_lattice(.01, 3);
    }
    if (with_non_corner_points) {
        params.set_quad_rule({
            { { { 1., 0., 0. } }, 0.1 },
            { { { 0., 1., 0. } }, 0.1 },
            { { { 0., 0., 1. } }, 0.1 },
            { { { .5, .5, 0. } }, 0.15 },
            { { { 0., .5, .5 } }, 0.15 },
            { { { .5, 0., .5 } }, 0.15 },
            { { { 1. / 3, 1. / 3, 1. / 3 } }, 0.25 },
        });
    }
    CAPTURE(barrier, params.area_weights, rule_kind);
    check_derivatives(s.mesh, s.X, params);
}

TEST_CASE(
    "ESP local transfer multiple EE samples and interior barrier",
    "[esp_potential][esp_local_transfer]")
{
    Eigen::MatrixXd V(6, 3);
    // Avoid endpoint feature switches: the potential is C1 there, so a unique
    // Hessian is not expected at the symmetric configuration.
    V << -.01, -.0003, .0006, .01, .0002, .0008, -.005, -.021, 0, -.0047, .019,
        0, .005, -.018, 0, .0053, .023, 0;
    Eigen::MatrixXi E(3, 2);
    E << 0, 1, 2, 3, 4, 5;
    CollisionMesh mesh(V, E);
    Eigen::VectorXi patches(6);
    patches << 0, 0, 1, 1, 2, 2;
    mesh.can_collide = make_vertex_patches_filter(patches);
    mesh.init_adjacencies();
    ESPParameters params(.03, .05);
    params.local_quadrature_transfer = true;
    check_derivatives(mesh, V, params);
    const ESPPotential p(params);
    double previous = -1;
    for (const double gap : { 1e-4, 1e-5, 1e-6, 1e-7 }) {
        Eigen::MatrixXd X = V;
        X.topRows(2).col(2).array() = gap;
        ESPCollisions c;
        c.build(mesh, X, params);
        const double energy = p(c, mesh, X);
        CHECK(std::isfinite(energy));
        CHECK(energy > previous);
        previous = energy;
    }
}

TEST_CASE(
    "ESP local transfer zero EE support and empty contacts",
    "[esp_potential][esp_local_transfer]")
{
    auto s = edge_over_box();
    ESPParameters legacy(.03, 0), local = legacy;
    local.local_quadrature_transfer = true;
    const ESPPotential a(legacy), b(local);
    ESPCollisions c;
    c.build(s.mesh, s.X, local);
    CHECK(std::abs(a(c, s.mesh, s.X) - b(c, s.mesh, s.X)) < 1e-14);
    CHECK(
        (a.gradient(c, s.mesh, s.X) - b.gradient(c, s.mesh, s.X)).norm()
        < 1e-12);
    ESPParameters params(.03, .05);
    params.local_quadrature_transfer = true;
    ESPPotential p(params);
    ESPCollisions empty;
    CHECK(p(empty, s.mesh, s.X) == 0);
    CHECK(p.gradient(empty, s.mesh, s.X).norm() == 0);
    CHECK(p.hessian(empty, s.mesh, s.X).nonZeros() == 0);
}

TEST_CASE(
    "ESP local transfer rejects incompatible quadrature modes",
    "[esp_potential][esp_local_transfer]")
{
    ESPParameters params(.03, .05);
    params.local_quadrature_transfer = true;
    CHECK_THROWS_AS(ESPPotential(params, false), std::invalid_argument);
    params.set_quad_rule({ { { { 1. / 3, 1. / 3, 1. / 3 } }, 1.0 } });
    CHECK_THROWS_AS(ESPPotential(params), std::invalid_argument);
}

TEST_CASE(
    "ESP local transfer rejects inconsistent friction weights",
    "[esp_potential][esp_local_transfer]")
{
    auto s = edge_over_box();
    ESPParameters params(.03, .05);
    params.local_quadrature_transfer = true;
    ESPCollisions c;
    c.build(s.mesh, s.X, params);
    TangentialCollisions friction;
    const Eigen::VectorXd zero = Eigen::VectorXd::Zero(s.X.rows());
    const Eigen::VectorXd mu = Eigen::VectorXd::Constant(s.X.rows(), .3);
    CHECK_NOTHROW(friction.build(s.mesh, s.X, c, params, 1.0, zero, zero));
    CHECK(friction.empty());
    CHECK_THROWS_AS(
        friction.build(s.mesh, s.X, c, params, 1.0, mu, mu),
        std::invalid_argument);
}
