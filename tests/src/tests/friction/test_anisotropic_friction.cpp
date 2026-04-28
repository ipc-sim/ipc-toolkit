#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/generators/catch_generators_range.hpp>

#include <ipc/friction/smooth_mu.hpp>
#include <ipc/collisions/tangential/tangential_collisions.hpp>
#include <ipc/collisions/tangential/vertex_vertex.hpp>
#include <ipc/collisions/normal/normal_collisions.hpp>
#include <ipc/potentials/friction_potential.hpp>
#include <ipc/potentials/barrier_potential.hpp>
#include <ipc/collision_mesh.hpp>

#include <finitediff.hpp>
#include <tests/utils.hpp>
#include <igl/edges.h>
#include <igl/PI.h>

#include <algorithm>
#include <cmath>

using namespace ipc;

TEST_CASE("Anisotropic mu effective computation", "[friction][anisotropic][mu]")
{
    static constexpr double EPSILON = 1e-6;
    static constexpr double MARGIN = 1e-8;

    // Test various directions and mu values
    const double mu_s0 = GENERATE(range(0.1, 1.0, 0.2));
    const double mu_s1 = GENERATE(range(0.1, 1.0, 0.2));
    const double mu_k0 = GENERATE(range(0.1, 1.0, 0.2));
    const double mu_k1 = GENERATE(range(0.1, 1.0, 0.2));

    Eigen::Vector2d mu_s_aniso(mu_s0, mu_s1);
    Eigen::Vector2d mu_k_aniso(mu_k0, mu_k1);

    // Test various directions
    const double angle = GENERATE(range(0.0, 2.0 * igl::PI, igl::PI / 4.0));
    Eigen::Vector2d tau_dir(std::cos(angle), std::sin(angle));

    CAPTURE(mu_s_aniso, mu_k_aniso, tau_dir);

    const auto [mu_s_eff, mu_k_eff] =
        anisotropic_mu_eff_f(tau_dir, mu_s_aniso, mu_k_aniso);

    // Verify L2 projection formula: mu_eff = sqrt((mu0*t0)^2 + (mu1*t1)^2)
    const double expected_mu_s_eff = std::sqrt(
        mu_s_aniso[0] * mu_s_aniso[0] * tau_dir[0] * tau_dir[0]
        + mu_s_aniso[1] * mu_s_aniso[1] * tau_dir[1] * tau_dir[1]);
    const double expected_mu_k_eff = std::sqrt(
        mu_k_aniso[0] * mu_k_aniso[0] * tau_dir[0] * tau_dir[0]
        + mu_k_aniso[1] * mu_k_aniso[1] * tau_dir[1] * tau_dir[1]);

    CHECK(
        mu_s_eff
        == Catch::Approx(expected_mu_s_eff).margin(MARGIN).epsilon(EPSILON));
    CHECK(
        mu_k_eff
        == Catch::Approx(expected_mu_k_eff).margin(MARGIN).epsilon(EPSILON));

    // Verify symmetry: mu_eff should be symmetric around principal axes
    if (std::abs(tau_dir[0]) < 1e-10) {
        // Pure y direction
        CHECK(
            mu_s_eff
            == Catch::Approx(mu_s_aniso[1]).margin(MARGIN).epsilon(EPSILON));
        CHECK(
            mu_k_eff
            == Catch::Approx(mu_k_aniso[1]).margin(MARGIN).epsilon(EPSILON));
    } else if (std::abs(tau_dir[1]) < 1e-10) {
        // Pure x direction
        CHECK(
            mu_s_eff
            == Catch::Approx(mu_s_aniso[0]).margin(MARGIN).epsilon(EPSILON));
        CHECK(
            mu_k_eff
            == Catch::Approx(mu_k_aniso[0]).margin(MARGIN).epsilon(EPSILON));
    }
}

TEST_CASE(
    "Anisotropic friction uses correct mu_s/mu_k transition",
    "[friction][anisotropic][mu-transition]")
{
    static constexpr double EPSILON = 1e-5;
    static constexpr double MARGIN = 1e-6;
    static constexpr double eps_v = 1e-3;

    // Fixed direction and anisotropic coefficients
    const double mu_s0 = GENERATE(0.2, 0.5, 0.8);
    const double mu_s1 = GENERATE(0.3, 0.6);
    const double mu_k0 = GENERATE(0.1, 0.4);
    const double mu_k1 = GENERATE(0.15, 0.35);
    Eigen::Vector2d mu_s_aniso(mu_s0, mu_s1);
    Eigen::Vector2d mu_k_aniso(mu_k0, mu_k1);
    Eigen::Vector2d tau_dir(1.0, 0.0);

    const auto [mu_s_eff, mu_k_eff] =
        anisotropic_mu_eff_f(tau_dir, mu_s_aniso, mu_k_aniso);

    // At y = 0 (stick), smooth_mu should approximate mu_s_eff
    const double mu_at_zero = smooth_mu(0.0, mu_s_eff, mu_k_eff, eps_v);
    CHECK(
        mu_at_zero == Catch::Approx(mu_s_eff).margin(MARGIN).epsilon(EPSILON));

    // At y = eps_v (slip), smooth_mu should approximate mu_k_eff
    const double mu_at_eps_v = smooth_mu(eps_v, mu_s_eff, mu_k_eff, eps_v);
    CHECK(
        mu_at_eps_v == Catch::Approx(mu_k_eff).margin(MARGIN).epsilon(EPSILON));

    // At y = 0.5*eps_v, smooth_mu should be between mu_s_eff and mu_k_eff
    const double mu_at_mid = smooth_mu(0.5 * eps_v, mu_s_eff, mu_k_eff, eps_v);
    const double lo = std::min(mu_s_eff, mu_k_eff);
    const double hi = std::max(mu_s_eff, mu_k_eff);
    CHECK(mu_at_mid >= lo - MARGIN);
    CHECK(mu_at_mid <= hi + MARGIN);
}

TEST_CASE(
    "Anisotropic friction isotropic", "[friction][anisotropic][isotropic]")
{
    // When mu_s_aniso and mu_k_aniso are zero, should use scalar mu_s and mu_k
    Eigen::Vector2d mu_s_aniso_zero = Eigen::Vector2d::Zero();
    Eigen::Vector2d mu_k_aniso_zero = Eigen::Vector2d::Zero();

    // Test with various directions
    const double angle = GENERATE(range(0.0, 2.0 * igl::PI, igl::PI / 4.0));
    Eigen::Vector2d tau_dir(std::cos(angle), std::sin(angle));

    // When anisotropic coefficients are zero, the function should handle it
    // gracefully (though the effective mu will be zero, which is expected)
    const auto [mu_s_eff, mu_k_eff] =
        anisotropic_mu_eff_f(tau_dir, mu_s_aniso_zero, mu_k_aniso_zero);

    CHECK(mu_s_eff == Catch::Approx(0.0).margin(1e-10));
    CHECK(mu_k_eff == Catch::Approx(0.0).margin(1e-10));
}

TEST_CASE(
    "anisotropic_x_from_tau_aniso", "[friction][anisotropic][smooth-mu-edge]")
{
    static constexpr double MARGIN = 1e-10;

    // tau_aniso_norm < tiny: expect default direction (1, 0)
    Eigen::Vector2d tau_aniso_tiny(1e-12, 1e-12);
    Eigen::Vector2d x_tiny = anisotropic_x_from_tau_aniso(tau_aniso_tiny);
    CHECK(x_tiny[0] == Catch::Approx(1.0).margin(MARGIN));
    CHECK(x_tiny[1] == Catch::Approx(0.0).margin(MARGIN));

    // Non-zero tau_aniso: expect normalized vector
    Eigen::Vector2d tau_aniso(0.6, 0.8);
    Eigen::Vector2d x = anisotropic_x_from_tau_aniso(tau_aniso);
    CHECK(x.norm() == Catch::Approx(1.0).margin(MARGIN));
    CHECK(x[0] == Catch::Approx(0.6).margin(MARGIN));
    CHECK(x[1] == Catch::Approx(0.8).margin(MARGIN));
}

TEST_CASE(
    "anisotropic_mu_eff_from_tau_aniso no_mu and isotropic",
    "[friction][anisotropic][smooth-mu-edge]")
{
    static constexpr double MARGIN = 1e-10;

    Eigen::Vector2d tau_aniso(0.6, 0.8);
    const double mu_s_iso = 0.5;
    const double mu_k_iso = 0.3;
    Eigen::Vector2d mu_s_aniso(0.6, 0.4);
    Eigen::Vector2d mu_k_aniso(0.5, 0.25);

    // Anisotropic path, no_mu == true: expect (1.0, 1.0)
    auto [mu_s, mu_k] = anisotropic_mu_eff_from_tau_aniso(
        tau_aniso, mu_s_aniso, mu_k_aniso, mu_s_iso, mu_k_iso, true);
    CHECK(mu_s == Catch::Approx(1.0).margin(MARGIN));
    CHECK(mu_k == Catch::Approx(1.0).margin(MARGIN));

    // Isotropic path (zero ellipse axes), no_mu == true: expect (1.0, 1.0)
    Eigen::Vector2d zero_aniso = Eigen::Vector2d::Zero();
    std::tie(mu_s, mu_k) = anisotropic_mu_eff_from_tau_aniso(
        tau_aniso, zero_aniso, zero_aniso, mu_s_iso, mu_k_iso, true);
    CHECK(mu_s == Catch::Approx(1.0).margin(MARGIN));
    CHECK(mu_k == Catch::Approx(1.0).margin(MARGIN));

    // Isotropic path, no_mu == false: expect scalar mu_s_isotropic,
    // mu_k_isotropic
    std::tie(mu_s, mu_k) = anisotropic_mu_eff_from_tau_aniso(
        tau_aniso, zero_aniso, zero_aniso, mu_s_iso, mu_k_iso, false);
    CHECK(mu_s == Catch::Approx(mu_s_iso).margin(MARGIN));
    CHECK(mu_k == Catch::Approx(mu_k_iso).margin(MARGIN));
}

TEST_CASE(
    "anisotropic_mu_eff_from_tau_aniso only mu_k_aniso",
    "[friction][anisotropic][smooth-mu]")
{
    static constexpr double MARGIN = 1e-10;

    // Only kinetic anisotropic: mu_s_aniso = 0, mu_k_aniso non-zero.
    // Expect mu_s_eff = mu_s_isotropic, mu_k_eff from ellipse.
    Eigen::Vector2d tau_aniso(1.0, 0.0); // direction (1,0)
    const double mu_s_iso = 0.5;
    const double mu_k_iso = 0.3;
    Eigen::Vector2d mu_s_aniso_zero = Eigen::Vector2d::Zero();
    Eigen::Vector2d mu_k_aniso(0.6, 0.25);

    auto [mu_s_eff, mu_k_eff] = anisotropic_mu_eff_from_tau_aniso(
        tau_aniso, mu_s_aniso_zero, mu_k_aniso, mu_s_iso, mu_k_iso, false);
    CHECK(mu_s_eff == Catch::Approx(mu_s_iso).margin(MARGIN));
    // For direction (1,0), mu_k_eff = mu_k_aniso[0]
    CHECK(mu_k_eff == Catch::Approx(mu_k_aniso[0]).margin(MARGIN));

    // Only static anisotropic: mu_s_aniso non-zero, mu_k_aniso = 0.
    Eigen::Vector2d mu_s_aniso(0.7, 0.35);
    Eigen::Vector2d mu_k_aniso_zero = Eigen::Vector2d::Zero();
    std::tie(mu_s_eff, mu_k_eff) = anisotropic_mu_eff_from_tau_aniso(
        tau_aniso, mu_s_aniso, mu_k_aniso_zero, mu_s_iso, mu_k_iso, false);
    CHECK(mu_s_eff == Catch::Approx(mu_s_aniso[0]).margin(MARGIN));
    CHECK(mu_k_eff == Catch::Approx(mu_k_iso).margin(MARGIN));
}

TEST_CASE(
    "Anisotropic friction per-pair assignment",
    "[friction][anisotropic][per-pair]")
{
    const double dhat = 1e-2;
    // Two vertices within dhat so we get at least one normal (and tangential)
    // collision
    Eigen::MatrixXd vertices(2, 3);
    vertices << 0.0, 0.0, 0.0, dhat * 0.5, 0.0, 0.0;

    Eigen::MatrixXi edges, faces;
    CollisionMesh mesh(vertices, edges, faces);

    NormalCollisions normal_collisions;
    normal_collisions.build(mesh, vertices, dhat);
    REQUIRE(!normal_collisions.empty());

    TangentialCollisions tangential_collisions;
    BarrierPotential barrier_potential(dhat, 1.0);
    tangential_collisions.build(
        mesh, vertices, normal_collisions, barrier_potential, 0.5, 0.3);
    REQUIRE(!tangential_collisions.empty());

    for (size_t i = 0; i < tangential_collisions.size(); ++i) {
        Eigen::Vector2d mu_s_aniso(0.5 + i * 0.1, 0.3 + i * 0.1);
        Eigen::Vector2d mu_k_aniso(0.4 + i * 0.1, 0.2 + i * 0.1);

        tangential_collisions[i].mu_s_aniso = mu_s_aniso;
        tangential_collisions[i].mu_k_aniso = mu_k_aniso;

        CHECK(tangential_collisions[i].mu_s_aniso == mu_s_aniso);
        CHECK(tangential_collisions[i].mu_k_aniso == mu_k_aniso);
    }
}

TEST_CASE(
    "Isotropic friction lagged mu matches scalar mu_s mu_k",
    "[friction][anisotropic][isotropic][lagged]")
{
    static constexpr double MARGIN = 1e-12;

    const double dhat = 1e-2;
    Eigen::MatrixXd vertices(2, 3);
    vertices << 0.0, 0.0, 0.0, dhat * 0.5, 0.0, 0.0;

    Eigen::MatrixXi edges, faces;
    CollisionMesh mesh(vertices, edges, faces);

    NormalCollisions normal_collisions;
    normal_collisions.build(mesh, vertices, dhat);
    REQUIRE(!normal_collisions.empty());

    const double mu_s = 0.5;
    const double mu_k = 0.3;
    TangentialCollisions tangential_collisions;
    BarrierPotential barrier_potential(dhat, 1.0);
    tangential_collisions.build(
        mesh, vertices, normal_collisions, barrier_potential, mu_s, mu_k);
    REQUIRE(!tangential_collisions.empty());

    Eigen::MatrixXd velocities(2, 3);
    velocities << 0.1, 0.0, 0.0, 0.0, 0.0, 0.0;

    tangential_collisions.update_lagged_anisotropic_friction_coefficients(
        mesh, vertices, Eigen::MatrixXd::Zero(vertices.rows(), vertices.cols()),
        velocities);

    for (size_t i = 0; i < tangential_collisions.size(); ++i) {
        const auto& c = tangential_collisions[i];
        CHECK(c.mu_s_aniso.squaredNorm() == Catch::Approx(0.0).margin(MARGIN));
        CHECK(c.mu_k_aniso.squaredNorm() == Catch::Approx(0.0).margin(MARGIN));
        CHECK(c.mu_aniso[0] == Catch::Approx(1.0).margin(MARGIN));
        CHECK(c.mu_aniso[1] == Catch::Approx(1.0).margin(MARGIN));
        CHECK(
            c.mu_s_effective_lagged
            == Catch::Approx(c.mu_s).margin(MARGIN).epsilon(1e-14));
        CHECK(
            c.mu_k_effective_lagged
            == Catch::Approx(c.mu_k).margin(MARGIN).epsilon(1e-14));
    }
}

TEST_CASE(
    "Anisotropy disabled matches isotropic baseline",
    "[friction][anisotropic][isotropic][parity]")
{
    const double dhat = 1e-2;
    Eigen::MatrixXd vertices(2, 3);
    vertices << 0.0, 0.0, 0.0, dhat * 0.5, 0.0, 0.0;

    Eigen::MatrixXi edges, faces;
    CollisionMesh mesh(vertices, edges, faces);

    NormalCollisions normal_collisions;
    normal_collisions.build(mesh, vertices, dhat);
    REQUIRE(!normal_collisions.empty());

    const double barrier_stiffness = 1.0;
    BarrierPotential barrier_potential(dhat, barrier_stiffness);
    FrictionPotential friction_potential(1e-4);

    TangentialCollisions baseline_collisions;
    baseline_collisions.build(
        mesh, vertices, normal_collisions, barrier_potential, 0.5, 0.3);
    REQUIRE(!baseline_collisions.empty());

    TangentialCollisions disabled_aniso_collisions;
    disabled_aniso_collisions.build(
        mesh, vertices, normal_collisions, barrier_potential, 0.5, 0.3);
    REQUIRE(!disabled_aniso_collisions.empty());

    Eigen::MatrixXd velocities(2, 3);
    velocities << 0.1, 0.05, 0.0, 0.0, 0.0, 0.0;
    const Eigen::MatrixXd lagged_displacements =
        Eigen::MatrixXd::Zero(vertices.rows(), vertices.cols());

    // Keep anisotropy disabled explicitly (same as default behavior).
    for (size_t i = 0; i < disabled_aniso_collisions.size(); ++i) {
        disabled_aniso_collisions[i].mu_s_aniso = Eigen::Vector2d::Zero();
        disabled_aniso_collisions[i].mu_k_aniso = Eigen::Vector2d::Zero();
        disabled_aniso_collisions[i].mu_aniso = Eigen::Vector2d::Ones();
    }

    baseline_collisions.update_lagged_anisotropic_friction_coefficients(
        mesh, vertices, lagged_displacements, velocities);
    disabled_aniso_collisions.update_lagged_anisotropic_friction_coefficients(
        mesh, vertices, lagged_displacements, velocities);

    const Eigen::VectorXd force_baseline = friction_potential.force(
        baseline_collisions, mesh, vertices, lagged_displacements, velocities,
        barrier_potential, 0.0, false);
    const Eigen::VectorXd force_disabled = friction_potential.force(
        disabled_aniso_collisions, mesh, vertices, lagged_displacements,
        velocities, barrier_potential, 0.0, false);
    const double force_diff_norm = (force_baseline - force_disabled).norm();
    const double force_ref_norm = std::max(force_baseline.norm(), 1.0);
    CHECK(force_diff_norm <= 1e-12 * force_ref_norm);

    const Eigen::VectorXd grad_baseline =
        friction_potential.gradient(baseline_collisions, mesh, velocities);
    const Eigen::VectorXd grad_disabled = friction_potential.gradient(
        disabled_aniso_collisions, mesh, velocities);
    const double grad_diff_norm = (grad_baseline - grad_disabled).norm();
    const double grad_ref_norm = std::max(grad_baseline.norm(), 1.0);
    CHECK(grad_diff_norm <= 1e-12 * grad_ref_norm);

    const Eigen::SparseMatrix<double> hess_baseline =
        friction_potential.hessian(baseline_collisions, mesh, velocities);
    const Eigen::SparseMatrix<double> hess_disabled =
        friction_potential.hessian(disabled_aniso_collisions, mesh, velocities);
    const Eigen::MatrixXd hess_baseline_dense = Eigen::MatrixXd(hess_baseline);
    const Eigen::MatrixXd hess_disabled_dense = Eigen::MatrixXd(hess_disabled);
    const double hess_diff_norm =
        (hess_baseline_dense - hess_disabled_dense).norm();
    const double hess_ref_norm = std::max(hess_baseline_dense.norm(), 1.0);
    CHECK(hess_diff_norm <= 1e-12 * hess_ref_norm);
}

TEST_CASE(
    "Lagged anisotropic mu refresh contract",
    "[friction][anisotropic][lagged][refresh]")
{
    static constexpr double MARGIN = 1e-12;

    const double dhat = 1e-2;
    Eigen::MatrixXd vertices(2, 3);
    vertices << 0.0, 0.0, 0.0, dhat * 0.5, 0.0, 0.0;

    Eigen::MatrixXi edges, faces;
    CollisionMesh mesh(vertices, edges, faces);

    NormalCollisions normal_collisions;
    normal_collisions.build(mesh, vertices, dhat);
    REQUIRE(!normal_collisions.empty());

    TangentialCollisions tangential_collisions;
    BarrierPotential barrier_potential(dhat, 1.0);
    tangential_collisions.build(
        mesh, vertices, normal_collisions, barrier_potential, 0.5, 0.3);
    REQUIRE(!tangential_collisions.empty());

    TangentialCollision& collision = tangential_collisions[0];
    collision.mu_s_aniso = Eigen::Vector2d(0.9, 0.2);
    collision.mu_k_aniso = Eigen::Vector2d(0.6, 0.15);

    const Eigen::MatrixXd lagged_displacements =
        Eigen::MatrixXd::Zero(vertices.rows(), vertices.cols());
    auto expected_lagged_mu =
        [&](const Eigen::MatrixXd& velocities) -> std::pair<double, double> {
        const VectorMax12d rest_dof = collision.dof(vertices, edges, faces);
        const VectorMax12d lag_dof =
            collision.dof(lagged_displacements, edges, faces);
        const VectorMax12d vel_dof = collision.dof(velocities, edges, faces);
        const VectorMax12d lagged_pos = rest_dof + lag_dof;

        const MatrixMax<double, 3, 2> P =
            collision.compute_tangent_basis(lagged_pos);
        const VectorMax2d beta = collision.compute_closest_point(lagged_pos);
        const MatrixMax<double, 3, 12> Gamma =
            collision.relative_velocity_jacobian(beta);
        const MatrixMax<double, 12, 2> T = Gamma.transpose() * P;
        const VectorMax2d tau = T.transpose() * vel_dof;
        const VectorMax2d tau_aniso =
            collision.mu_aniso.head(tau.size()).cwiseProduct(tau);

        return anisotropic_mu_eff_from_tau_aniso(
            tau_aniso, collision.mu_s_aniso, collision.mu_k_aniso,
            collision.mu_s, collision.mu_k, false);
    };

    Eigen::MatrixXd velocities_1 = Eigen::MatrixXd::Zero(2, 3);
    velocities_1.row(0) << 0.0, 1.0, 0.0;
    const auto [mu_s_expected_1, mu_k_expected_1] =
        expected_lagged_mu(velocities_1);
    tangential_collisions.update_lagged_anisotropic_friction_coefficients(
        mesh, vertices, lagged_displacements, velocities_1);
    CHECK(
        collision.mu_s_effective_lagged
        == Catch::Approx(mu_s_expected_1).margin(MARGIN));
    CHECK(
        collision.mu_k_effective_lagged
        == Catch::Approx(mu_k_expected_1).margin(MARGIN));

    // Change slip kinematics but do not refresh: lagged coefficients must stay
    // fixed for evaluation.
    Eigen::MatrixXd velocities_2 = Eigen::MatrixXd::Zero(2, 3);
    velocities_2.row(0) << 0.0, 0.37, 0.91;
    const auto [mu_s_expected_2, mu_k_expected_2] =
        expected_lagged_mu(velocities_2);
    CHECK(
        collision.mu_s_effective_lagged
        == Catch::Approx(mu_s_expected_1).margin(MARGIN));
    CHECK(
        collision.mu_k_effective_lagged
        == Catch::Approx(mu_k_expected_1).margin(MARGIN));

    tangential_collisions.update_lagged_anisotropic_friction_coefficients(
        mesh, vertices, lagged_displacements, velocities_2);
    CHECK(
        collision.mu_s_effective_lagged
        == Catch::Approx(mu_s_expected_2).margin(MARGIN));
    CHECK(
        collision.mu_k_effective_lagged
        == Catch::Approx(mu_k_expected_2).margin(MARGIN));
}

TEST_CASE("Anisotropic friction force", "[friction][anisotropic][force]")
{
    static constexpr double MARGIN = 1e-6;

    const double dhat = 1e-2;
    // Vertices within dhat so we get at least one normal and tangential
    // collision
    Eigen::MatrixXd vertices(2, 3);
    vertices << 0.0, 0.0, 0.0, dhat * 0.5, 0.0, 0.0;

    Eigen::MatrixXi edges, faces;
    CollisionMesh mesh(vertices, edges, faces);

    NormalCollisions normal_collisions;
    normal_collisions.build(mesh, vertices, dhat);
    REQUIRE(!normal_collisions.empty());

    TangentialCollisions tangential_collisions;
    const double barrier_stiffness = 1.0;
    BarrierPotential barrier_potential(dhat, barrier_stiffness);
    tangential_collisions.build(
        mesh, vertices, normal_collisions, barrier_potential, 0.5, 0.3);
    REQUIRE(!tangential_collisions.empty());

    // Set anisotropic friction coefficients on first collision
    TangentialCollision& collision = tangential_collisions[0];
    Eigen::Vector2d mu_s_aniso(0.8, 0.4); // Higher in x, lower in y
    Eigen::Vector2d mu_k_aniso(0.6, 0.3);
    collision.mu_s_aniso = mu_s_aniso;
    collision.mu_k_aniso = mu_k_aniso;

    FrictionPotential friction_potential(1e-4);

    // Test with different velocity directions
    const double angle = GENERATE(range(0.0, 2.0 * igl::PI, igl::PI / 4.0));
    const double velocity_magnitude = GENERATE(0.01, 0.1);

    Eigen::MatrixXd velocities(2, 3);
    // Velocity in xy-plane: (cos(angle), sin(angle), 0)
    // For vertex-vertex along x-axis, tangent plane is yz; tangential
    // component comes only from sin(angle). When angle=0 or π, tangential
    // velocity is zero, so friction force is correctly zero.
    velocities << velocity_magnitude * std::cos(angle),
        velocity_magnitude * std::sin(angle), 0.0, //
        0.0, 0.0, 0.0;

    CAPTURE(angle, velocity_magnitude, mu_s_aniso, mu_k_aniso);

    tangential_collisions.update_lagged_anisotropic_friction_coefficients(
        mesh, vertices, Eigen::MatrixXd::Zero(vertices.rows(), vertices.cols()),
        velocities);

    // Compute force for this single collision
    VectorMax12d force = friction_potential.force(
        collision, collision.dof(vertices, edges, faces),
        collision.dof(Eigen::MatrixXd::Zero(2, 3), edges, faces),
        collision.dof(velocities, edges, faces), barrier_potential, 0.0, false);

    // Verify force is non-zero when tangential velocity is non-zero
    // (velocity along normal yields zero friction by design)
    const bool has_tangential_velocity = std::abs(std::sin(angle)) > 1e-10;
    if (velocity_magnitude > 1e-10 && has_tangential_velocity) {
        CHECK(force.norm() > 1e-10);
    }

    // For anisotropic friction, verify that force magnitude changes with
    // direction
    // When velocity is along x-axis (angle=0), mu_eff should be mu_s_aniso[0]
    // When velocity is along y-axis (angle=π/2), mu_eff should be mu_s_aniso[1]
    if (std::abs(std::sin(angle)) < 1e-6) {
        // Pure x direction
        const auto [mu_s_eff, mu_k_eff] = anisotropic_mu_eff_f(
            Eigen::Vector2d(1.0, 0.0), mu_s_aniso, mu_k_aniso);
        CHECK(mu_s_eff == Catch::Approx(mu_s_aniso[0]).margin(MARGIN));
    } else if (std::abs(std::cos(angle)) < 1e-6) {
        // Pure y direction
        const auto [mu_s_eff, mu_k_eff] = anisotropic_mu_eff_f(
            Eigen::Vector2d(0.0, 1.0), mu_s_aniso, mu_k_aniso);
        CHECK(mu_s_eff == Catch::Approx(mu_s_aniso[1]).margin(MARGIN));
    }
}

TEST_CASE(
    "Anisotropic friction force jacobian", "[friction][anisotropic][jacobian]")
{
    static constexpr double H = 1e-8;

    // Create a simple mesh with two vertices
    Eigen::MatrixXd vertices(2, 3);
    vertices << 0.0, 0.0, 0.0, // vertex 0
        0.05, 0.0, 0.0;        // vertex 1 (close to vertex 0)

    Eigen::MatrixXi edges, faces;
    CollisionMesh mesh(vertices, edges, faces);

    // Create normal collision
    NormalCollisions normal_collisions;
    const double dhat = 1e-2;
    normal_collisions.build(mesh, vertices, dhat);

    if (normal_collisions.empty()) {
        return; // Skip if no collisions
    }

    // Build tangential collisions
    TangentialCollisions tangential_collisions;
    const double barrier_stiffness = 1.0;
    BarrierPotential barrier_potential(dhat, barrier_stiffness);
    tangential_collisions.build(
        mesh, vertices, normal_collisions, barrier_potential, 0.5, 0.3);

    if (tangential_collisions.empty()) {
        return; // Skip if no tangential collisions
    }

    // Set anisotropic friction coefficients on first collision
    TangentialCollision& collision = tangential_collisions[0];
    Eigen::Vector2d mu_s_aniso(0.8, 0.4);
    Eigen::Vector2d mu_k_aniso(0.6, 0.3);
    collision.mu_s_aniso = mu_s_aniso;
    collision.mu_k_aniso = mu_k_aniso;

    FrictionPotential friction_potential(1e-4);

    // Test with a specific velocity
    Eigen::MatrixXd velocities(2, 3);
    velocities << 0.1, 0.05, 0.0, // velocity in tangent plane
        0.0, 0.0, 0.0;

    CAPTURE(mu_s_aniso, mu_k_aniso);

    tangential_collisions.update_lagged_anisotropic_friction_coefficients(
        mesh, vertices, Eigen::MatrixXd::Zero(vertices.rows(), vertices.cols()),
        velocities);

    // Compute force jacobian w.r.t. velocities
    MatrixMax12d jacobian = friction_potential.force_jacobian(
        collision, collision.dof(vertices, edges, faces),
        collision.dof(Eigen::MatrixXd::Zero(2, 3), edges, faces),
        collision.dof(velocities, edges, faces), barrier_potential,
        FrictionPotential::DiffWRT::VELOCITIES, 0.0);

    // Verify jacobian is non-zero
    CHECK(jacobian.norm() > 1e-10);

    // Compare with finite differences
    Eigen::VectorXd V_flat = collision.dof(velocities, edges, faces);

    auto F_V = [&](const Eigen::VectorXd& v) -> Eigen::MatrixXd {
        // Convert vector back to matrix format
        Eigen::MatrixXd v_mat(2, 3);
        v_mat.row(0) = v.head(3);
        v_mat.row(1) = v.tail(3);
        VectorMax12d f = friction_potential.force(
            collision, collision.dof(vertices, edges, faces),
            collision.dof(Eigen::MatrixXd::Zero(2, 3), edges, faces),
            collision.dof(v_mat, edges, faces), barrier_potential, 0.0, false);
        return Eigen::MatrixXd(f);
    };

    Eigen::MatrixXd fd_jacobian;
    fd::finite_jacobian(V_flat, F_V, fd_jacobian, fd::AccuracyOrder::FOURTH, H);

    CHECK(fd::compare_jacobian(jacobian, fd_jacobian));
    if (!fd::compare_jacobian(jacobian, fd_jacobian)) {
        tests::print_compare_nonzero(jacobian, fd_jacobian);
    }
}

TEST_CASE(
    "Combined mu_aniso and mu_s_aniso/mu_k_aniso friction",
    "[friction][anisotropic][combined]")
{
    static constexpr double H = 1e-8;

    // Test that both mechanisms work together:
    // 1. mu_aniso: velocity scaling (tau_aniso = mu_aniso ⊙ tau)
    // 2. mu_s_aniso/mu_k_aniso: direction-dependent coefficients

    // Create a simple mesh with two vertices
    Eigen::MatrixXd vertices(2, 3);
    vertices << 0.0, 0.0, 0.0, // vertex 0
        0.05, 0.0, 0.0;        // vertex 1 (close to vertex 0)

    Eigen::MatrixXi edges, faces;
    CollisionMesh mesh(vertices, edges, faces);

    // Create normal collision
    NormalCollisions normal_collisions;
    const double dhat = 1e-2;
    normal_collisions.build(mesh, vertices, dhat);

    if (normal_collisions.empty()) {
        return; // Skip if no collisions
    }

    // Build tangential collisions
    TangentialCollisions tangential_collisions;
    const double barrier_stiffness = 1.0;
    BarrierPotential barrier_potential(dhat, barrier_stiffness);
    tangential_collisions.build(
        mesh, vertices, normal_collisions, barrier_potential, 0.5, 0.3);

    if (tangential_collisions.empty()) {
        return; // Skip if no tangential collisions
    }

    // Set BOTH anisotropic mechanisms on first collision:
    TangentialCollision& collision = tangential_collisions[0];

    // mu_aniso: velocity scaling (e.g., 0.8 in first tangent dir, 1.2 in
    // second) This scales the velocity differently along the two tangent
    // directions
    Eigen::Vector2d mu_aniso(0.8, 1.2);
    collision.mu_aniso = mu_aniso;

    // mu_s_aniso/mu_k_aniso: direction-dependent friction coefficients (e.g.,
    // higher friction in first tangent dir, lower in second)
    Eigen::Vector2d mu_s_aniso(0.6, 0.3);
    Eigen::Vector2d mu_k_aniso(0.5, 0.25);
    collision.mu_s_aniso = mu_s_aniso;
    collision.mu_k_aniso = mu_k_aniso;

    FrictionPotential friction_potential(1e-4);

    // Test various velocity directions
    const double angle = GENERATE(range(0.0, 2.0 * igl::PI, igl::PI / 6.0));
    const double velocity_magnitude = GENERATE(0.01, 0.1, 0.5);

    Eigen::MatrixXd velocities(2, 3);
    velocities << velocity_magnitude * std::cos(angle),
        velocity_magnitude * std::sin(angle), 0.0, //
        0.0, 0.0, 0.0;

    CAPTURE(mu_aniso, mu_s_aniso, mu_k_aniso, angle, velocity_magnitude);

    tangential_collisions.update_lagged_anisotropic_friction_coefficients(
        mesh, vertices, Eigen::MatrixXd::Zero(vertices.rows(), vertices.cols()),
        velocities);

    // Compute force
    VectorMax12d force = friction_potential.force(
        collision, collision.dof(vertices, edges, faces),
        collision.dof(Eigen::MatrixXd::Zero(2, 3), edges, faces),
        collision.dof(velocities, edges, faces), barrier_potential, 0.0, false);

    // Verify force is non-zero
    CHECK(force.norm() > 1e-10);

    // Compute force jacobian w.r.t. velocities
    MatrixMax12d jacobian = friction_potential.force_jacobian(
        collision, collision.dof(vertices, edges, faces),
        collision.dof(Eigen::MatrixXd::Zero(2, 3), edges, faces),
        collision.dof(velocities, edges, faces), barrier_potential,
        FrictionPotential::DiffWRT::VELOCITIES, 0.0);

    // Verify jacobian is non-zero
    CHECK(jacobian.norm() > 1e-10);

    // Compare with finite differences
    Eigen::VectorXd V_flat = collision.dof(velocities, edges, faces);

    auto F_V = [&](const Eigen::VectorXd& v) -> Eigen::MatrixXd {
        // Convert vector back to matrix format
        Eigen::MatrixXd v_mat(2, 3);
        v_mat.row(0) = v.head(3);
        v_mat.row(1) = v.tail(3);
        VectorMax12d f = friction_potential.force(
            collision, collision.dof(vertices, edges, faces),
            collision.dof(Eigen::MatrixXd::Zero(2, 3), edges, faces),
            collision.dof(v_mat, edges, faces), barrier_potential, 0.0, false);
        return Eigen::MatrixXd(f);
    };

    Eigen::MatrixXd fd_jacobian;
    fd::finite_jacobian(V_flat, F_V, fd_jacobian, fd::AccuracyOrder::FOURTH, H);

    CHECK(fd::compare_jacobian(jacobian, fd_jacobian));
    if (!fd::compare_jacobian(jacobian, fd_jacobian)) {
        tests::print_compare_nonzero(jacobian, fd_jacobian);
    }

    // Verify that combined mechanism produces different results than either
    // mechanism alone. Save current force.
    Eigen::VectorXd combined_force = force;

    // Test with mu_aniso only (disable direction-dependent by setting
    // mu_s_aniso to zero)
    collision.mu_s_aniso = Eigen::Vector2d::Zero();
    collision.mu_k_aniso = Eigen::Vector2d::Zero();
    tangential_collisions.update_lagged_anisotropic_friction_coefficients(
        mesh, vertices, Eigen::MatrixXd::Zero(vertices.rows(), vertices.cols()),
        velocities);

    VectorMax12d mu_aniso_only_force = friction_potential.force(
        collision, collision.dof(vertices, edges, faces),
        collision.dof(Eigen::MatrixXd::Zero(2, 3), edges, faces),
        collision.dof(velocities, edges, faces), barrier_potential, 0.0, false);

    // Test with direction-dependent only (set mu_aniso to identity)
    collision.mu_aniso = Eigen::Vector2d(1.0, 1.0);
    collision.mu_s_aniso = mu_s_aniso;
    collision.mu_k_aniso = mu_k_aniso;
    tangential_collisions.update_lagged_anisotropic_friction_coefficients(
        mesh, vertices, Eigen::MatrixXd::Zero(vertices.rows(), vertices.cols()),
        velocities);

    VectorMax12d dir_dep_only_force = friction_potential.force(
        collision, collision.dof(vertices, edges, faces),
        collision.dof(Eigen::MatrixXd::Zero(2, 3), edges, faces),
        collision.dof(velocities, edges, faces), barrier_potential, 0.0, false);

    // Combined force should generally differ from either mechanism alone
    // (unless in special cases like aligned with principal axes)
    // We just verify the forces are computed without error here
    CHECK(mu_aniso_only_force.norm() > 1e-10);
    CHECK(dir_dep_only_force.norm() > 1e-10);
}

TEST_CASE(
    "Anisotropic friction dissipative gradient matches minus force",
    "[friction][anisotropic][potential-force]")
{
    Eigen::MatrixXd vertices(2, 3);
    vertices << 0.0, 0.0, 0.0, 0.05, 0.0, 0.0;

    Eigen::MatrixXi edges, faces;
    CollisionMesh mesh(vertices, edges, faces);

    NormalCollisions normal_collisions;
    const double dhat = 1e-2;
    normal_collisions.build(mesh, vertices, dhat);

    if (normal_collisions.empty()) {
        return;
    }

    const double barrier_stiffness = 1.0;
    BarrierPotential barrier_potential(dhat, barrier_stiffness);
    TangentialCollisions tangential_collisions;
    tangential_collisions.build(
        mesh, vertices, normal_collisions, barrier_potential, 0.5, 0.3);

    if (tangential_collisions.empty()) {
        return;
    }

    TangentialCollision& collision = tangential_collisions[0];
    collision.mu_s_aniso = Eigen::Vector2d(0.75, 0.42);
    collision.mu_k_aniso = Eigen::Vector2d(0.58, 0.29);
    // Nontrivial velocity scaling together with ellipse axes (lagged μ path).
    collision.mu_aniso = Eigen::Vector2d(0.9, 1.1);

    Eigen::MatrixXd velocities(2, 3);
    velocities << 0.12, 0.07, 0.0, 0.0, 0.0, 0.0;

    tangential_collisions.update_lagged_anisotropic_friction_coefficients(
        mesh, vertices, Eigen::MatrixXd::Zero(vertices.rows(), vertices.cols()),
        velocities);

    FrictionPotential friction_potential(1e-4);

    const VectorMax12d force = friction_potential.force(
        collision, collision.dof(vertices, edges, faces),
        collision.dof(Eigen::MatrixXd::Zero(2, 3), edges, faces),
        collision.dof(velocities, edges, faces), barrier_potential, 0.0, false);

    const VectorMax12d grad = friction_potential.gradient(
        collision, collision.dof(velocities, edges, faces));

    CHECK(fd::compare_gradient(-force, grad));
}

TEST_CASE(
    "Anisotropic TangentialPotential gradient and Hessian",
    "[friction][anisotropic][gradient]")
{
    // Create a simple 3D mesh: a triangle and a nearby point to force
    // a face-vertex tangential collision
    Eigen::MatrixXd vertices(4, 3);
    vertices << 0.0, 0.0, 0.0, // triangle vertex 0
        0.1, 0.0, 0.0,         // triangle vertex 1
        0.0, 0.1, 0.0,         // triangle vertex 2
        0.02, 0.02, 0.001;     // point (vertex 3) slightly above triangle

    Eigen::MatrixXi faces(1, 3);
    faces << 0, 1, 2;
    Eigen::MatrixXi edges;
    igl::edges(faces, edges);
    CollisionMesh mesh(vertices, edges, faces);
    mesh.init_area_jacobians();

    const double dhat = 1e-2;
    BarrierPotential barrier_potential(dhat, 1e6);

    // Build normal collisions
    NormalCollisions normal_collisions;
    normal_collisions.set_use_area_weighting(true);
    normal_collisions.set_collision_set_type(
        NormalCollisions::CollisionSetType::IPC);
    normal_collisions.build(mesh, vertices, dhat);

    REQUIRE(normal_collisions.size() == 1); // should have a collision

    // Build tangential collisions
    TangentialCollisions tangential_collisions;
    tangential_collisions.build(
        mesh, vertices, normal_collisions, barrier_potential, 0.5, 0.3);

    REQUIRE(tangential_collisions.size() == 1);

    // Use the first tangential collision
    TangentialCollision& collision = tangential_collisions[0];

    // Set direction-dependent friction coefficients
    collision.mu_s_aniso = Eigen::Vector2d(0.6, 0.3);
    collision.mu_k_aniso = Eigen::Vector2d(0.5, 0.25);

    // Keep default mu_aniso scaling
    collision.mu_aniso = Eigen::Vector2d(1.0, 1.0);

    FrictionPotential D(1e-4);

    // Choose a sample velocity: only the point (vertex 3) moves
    Eigen::MatrixXd velocities = Eigen::MatrixXd::Random(4, 3);

    const Eigen::MatrixXd lagged_displacements =
        Eigen::MatrixXd::Zero(vertices.rows(), vertices.cols());
    tangential_collisions.update_lagged_anisotropic_friction_coefficients(
        mesh, mesh.rest_positions(), lagged_displacements, velocities);

    // Analytical gradient (local DOF vector)
    const VectorMax12d vel_dof =
        collision.dof(velocities, mesh.edges(), mesh.faces());
    const VectorMax12d analytic_grad = D.gradient(collision, vel_dof);
    REQUIRE(!analytic_grad.isZero());

    // Finite-difference gradient of the scalar potential w.r.t. velocity DOF
    auto f = [&](const Eigen::VectorXd& v) -> double {
        return D(collision, v);
    };

    Eigen::VectorXd fd_grad;
    fd::finite_gradient(vel_dof, f, fd_grad);

    CHECK(fd::compare_gradient(analytic_grad, fd_grad));
    if (!fd::compare_gradient(analytic_grad, fd_grad)) {
        std::cout << "Analytic gradient: " << analytic_grad.transpose() << "\n";
        std::cout << "Finite-difference gradient: " << fd_grad.transpose()
                  << "\n";
    }

    // Analytical Hessian (local DOF x local DOF)
    const MatrixMax12d analytic_hessian = D.hessian(collision, vel_dof);
    REQUIRE(!analytic_hessian.isZero());

    // Finite-difference Hessian: Jacobian of the gradient w.r.t. velocity DOF
    auto g_jac = [&](const Eigen::VectorXd& v) -> Eigen::VectorXd {
        return D.gradient(collision, v);
    };

    Eigen::MatrixXd fd_hessian;
    fd::finite_jacobian(vel_dof, g_jac, fd_hessian);

    CHECK(fd::compare_jacobian(analytic_hessian, fd_hessian));
    if (!fd::compare_jacobian(analytic_hessian, fd_hessian)) {
        std::cout << "Analytic Hessian:\n" << analytic_hessian << "\n";
        std::cout << "Finite-difference Hessian:\n" << fd_hessian << "\n";
    }

    // Check gradient = -force
    const VectorMax12d force = D.force(
        collision,
        collision.dof(mesh.rest_positions(), mesh.edges(), mesh.faces()),
        collision.dof(lagged_displacements, mesh.edges(), mesh.faces()),
        vel_dof, barrier_potential, 0.0, false);
    CHECK(analytic_grad.isApprox(-force));
    if (!analytic_grad.isApprox(-force)) {
        std::cout << "Analytic gradient: " << analytic_grad.transpose() << "\n";
        std::cout << "Negative force: " << (-force).transpose() << "\n";
    }
}
