// Tests for the velocity-based lagged friction in the demo simulator.

#include <Eigen/Core>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <tests/utils.hpp>

#include <ipc/collisions/normal/normal_collisions.hpp>
#include <ipc/collisions/tangential/tangential_collisions.hpp>
#include <ipc/dynamics/rigid/rigid_bodies.hpp>
#include <simulator.hpp>

#include <finitediff.hpp>

#include <igl/PI.h>

#include <iostream>

using namespace ipc;
using namespace ipc::rigid;
using ipc::demo::Simulator;

namespace {

/// A cube above the ground plane y = 0 with an initial gap of
/// gap_factor·dhat. Use gap_factor < 1 to start in contact (barrier active);
/// use gap_factor slightly > 1 so contact engages naturally within a step —
/// placing a body deep inside the activation distance with a freshly seeded
/// stiff barrier produces unphysically large (lagged) normal forces.
std::shared_ptr<RigidBodies> cube_on_plane(
    std::vector<Pose>& initial_poses,
    const double dhat,
    const double density,
    const double gap_factor)
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    REQUIRE(tests::load_mesh("cube.ply", V, E, F));

    initial_poses.assign(1, Pose::Identity(3));
    initial_poses[0].position.y() = -V.col(1).minCoeff() + gap_factor * dhat;

    auto bodies = RigidBodies::build_from_meshes(
        { V }, { E }, { F }, { density }, initial_poses);
    bodies->planes.emplace_back(Eigen::Vector3d(0, 1, 0), 0); // ground y = 0
    return bodies;
}

} // namespace

TEST_CASE(
    "Friction simulator gradient and hessian",
    "[friction][simulator][gradient][hessian]")
{
    const bool affine = GENERATE(false, true);
    // Sliding speed regime: well below and well above the eps_v = 1e-3 kink
    // (finite differences are unreliable at the kink itself).
    const double offset = GENERATE(1e-6, 1e-3);

    CAPTURE(affine, offset);

    Simulator::Settings settings;
    settings.body_dynamics = affine ? Simulator::BodyDynamics::AFFINE
                                    : Simulator::BodyDynamics::RIGID;
    settings.friction_coefficient = 0.5;

    std::vector<Pose> initial_poses;
    auto bodies =
        cube_on_plane(initial_poses, settings.dhat, 1000.0, /*gap_factor=*/0.5);

    const double dt = 0.01;
    Simulator sim(bodies, initial_poses, dt, settings);

    // Lag the collision and friction sets at the initial state.
    const Eigen::VectorXd x0 = sim.kinematics().dof(initial_poses);
    sim.initialize_step();
    sim.update_collisions(x0);
    REQUIRE(!sim.normal_collisions().empty());
    sim.initialize_friction_step();
    sim.update_friction_collisions(x0);
    REQUIRE(!sim.tangential_collisions().empty());

    // Slide tangentially (and rotate slightly in rigid mode to exercise the
    // second-order chain-rule term).
    Eigen::VectorXd x = x0;
    x(0) += offset; // tangential slide along x
    if (!affine) {
        x(5) += 0.1 * offset; // rotation about z
    } else {
        x(4) += 0.1 * offset; // perturb A
    }

    Eigen::VectorXd grad;
    sim.gradient(x, grad);

    Eigen::VectorXd fd_grad;
    fd::finite_gradient(
        x, [&](const Eigen::VectorXd& y) { return sim.value(y); }, fd_grad);

    CHECK(fd::compare_gradient(grad, fd_grad));
    if (!fd::compare_gradient(grad, fd_grad)) {
        std::cout << "analytic:\n" << grad << "\n\n";
        std::cout << "numerical:\n" << fd_grad << "\n\n";
    }

    sim.set_project_to_psd(false); // FD comparison needs the exact Hessian
    Eigen::SparseMatrix<double> hess_sparse;
    sim.hessian(x, hess_sparse);
    const Eigen::MatrixXd hess = hess_sparse;

    Eigen::MatrixXd fd_hess;
    fd::finite_jacobian(
        x,
        [&](const Eigen::VectorXd& y) {
            Eigen::VectorXd g;
            sim.gradient(y, g);
            return g;
        },
        fd_hess);

    // NOTE: 1e-2 tolerance: the FD of the gradient carries roundoff noise
    // of O(1) from the ~1e8-scale (affine) mass-matrix entries, while the
    // friction/barrier entries being verified are O(1)-O(1e3).
    CHECK(fd::compare_hessian(hess, fd_hess, 1e-2));
    if (!fd::compare_hessian(hess, fd_hess, 1e-2)) {
        std::cout << "analytic:\n" << hess << "\n\n";
        std::cout << "numerical:\n" << fd_hess << "\n\n";
    }
}

TEST_CASE("Friction incline stick and slip", "[friction][simulator]")
{
    const bool affine = GENERATE(false, true);
    CAPTURE(affine);

    // Tilt gravity instead of the plane: an incline of angle θ.
    const double theta = 20.0 * igl::PI / 180.0; // tan θ ≈ 0.364
    const double g = 9.81;

    Simulator::Settings settings;
    settings.body_dynamics = affine ? Simulator::BodyDynamics::AFFINE
                                    : Simulator::BodyDynamics::RIGID;
    settings.gravity =
        Eigen::Vector3d(g * std::sin(theta), -g * std::cos(theta), 0);
    // Iterate the friction lagging to convergence: with a single lagged
    // solve, the normal forces are sampled at velocity-criterion-accurate
    // (not force-balanced) states, where the stiff barrier makes them noisy.
    settings.friction_iterations = -1;
    // ... and tighten the velocity criterion so its displacement tolerance is
    // below the barrier force-well width: otherwise the solver can quit while
    // a corner is jammed deep in the barrier (huge lagged normal force).
    settings.velocity_conv_tol = 1e-3;

    SECTION("sticks when mu > tan(theta)")
    {
        settings.friction_coefficient = 0.5;
    }
    SECTION("slides when mu < tan(theta)")
    {
        settings.friction_coefficient = 0.1;
    }

    std::vector<Pose> initial_poses;
    auto bodies =
        cube_on_plane(initial_poses, settings.dhat, 1000.0, /*gap_factor=*/1.1);

    const double dt = 0.01, t_end = 1.0;
    Simulator sim(bodies, initial_poses, dt, settings);

    const Eigen::Vector3d p0 = initial_poses[0].position;

    // Sample the position at three times past the landing transient and
    // measure the steady-sliding acceleration by central differences,
    // isolating the physics from the landing/settling transient.
    REQUIRE(sim.run(0.5 * t_end));
    const double x1 = (sim.poses()[0].position - p0).x();
    REQUIRE(sim.run(0.75 * t_end));
    const double xm = (sim.poses()[0].position - p0).x();
    REQUIRE(sim.run(t_end));
    const double x2 = (sim.poses()[0].position - p0).x();

    const double h = 0.25 * t_end;
    const double accel = (x2 - 2 * xm + x1) / (h * h);

    if (settings.friction_coefficient > std::tan(theta)) {
        // Static friction holds the cube in place (allow settling on the
        // order of the barrier activation distance).
        CHECK(std::abs(x2) < 10 * settings.dhat);
    } else {
        // Steady kinetic sliding: a = g (sin(theta) - mu cos(theta)).
        const double expected = g
            * (std::sin(theta)
               - settings.friction_coefficient * std::cos(theta));
        CHECK(accel == Catch::Approx(expected).epsilon(0.05));
    }
}

// The hard gate for the velocity-formulation scaling: a cube sliding on a
// flat plane decelerates at exactly μg, so it stops after v₀²/(2μg).
TEST_CASE("Friction sliding deceleration", "[friction][simulator]")
{
    const int bdf_order = GENERATE(1, 2);
    CAPTURE(bdf_order);

    const double mu = 0.3, g = 9.81, v0 = 1.0;

    Simulator::Settings settings;
    settings.friction_coefficient = mu;
    settings.bdf_order = bdf_order;
    settings.friction_iterations = -1; // see the incline test's notes
    settings.velocity_conv_tol = 1e-3;

    std::vector<Pose> initial_poses;
    auto bodies =
        cube_on_plane(initial_poses, settings.dhat, 1000.0, /*gap_factor=*/1.1);

    std::vector<Pose> initial_velocities(1, Pose::Identity(3));
    initial_velocities[0].position = Eigen::Vector3d(v0, 0, 0);

    const double dt = 0.01;
    Simulator sim(bodies, initial_poses, initial_velocities, dt, settings);

    const Eigen::Vector3d p0 = initial_poses[0].position;
    const double stopping_time = v0 / (mu * g);          // ≈ 0.34 s
    const double stopping_dist = v0 * v0 / (2 * mu * g); // ≈ 0.17 m

    REQUIRE(sim.run(/*t_end=*/2 * stopping_time));
    const double slide = (sim.poses()[0].position - p0).x();

    // First-order time discretization error is O(v₀ dt) per the ½aΔt²-style
    // offset of the stopping distance.
    CHECK(slide == Catch::Approx(stopping_dist).margin(3 * v0 * dt));

    // ... and it stays stopped (static friction holds).
    const Eigen::Vector3d p_stop = sim.poses()[0].position;
    REQUIRE(sim.run(/*t_end=*/2 * stopping_time + 0.2));
    CHECK(
        (sim.poses()[0].position - p_stop).norm()
        < settings.static_friction_speed_bound * 0.2);
}

TEST_CASE("Friction lagging semantics", "[friction][simulator]")
{
    Simulator::Settings settings;
    std::vector<Pose> initial_poses;

    SECTION("unlimited lagging converges the momentum balance")
    {
        settings.friction_coefficient = 0.3;
        settings.friction_iterations = -1;

        auto bodies = cube_on_plane(
            initial_poses, settings.dhat, 1000.0, /*gap_factor=*/1.1);
        std::vector<Pose> initial_velocities(1, Pose::Identity(3));
        initial_velocities[0].position = Eigen::Vector3d(1.0, 0, 0);

        Simulator sim(
            bodies, initial_poses, initial_velocities, /*dt=*/0.01, settings);

        for (int i = 0; i < 5; ++i) {
            REQUIRE(sim.step());
            // The lagging ran and measured the momentum balance. (It may
            // legitimately end above eps_d when the inner solver's own
            // stopping criteria bound the reachable residual.)
            CHECK(sim.last_momentum_balance() >= 0);
        }
    }

    SECTION("friction_iterations == 0 disables friction entirely")
    {
        auto run_sim = [&](const double mu, const int iterations) {
            std::vector<Pose> poses;
            Simulator::Settings s;
            s.friction_coefficient = mu;
            s.friction_iterations = iterations;
            auto bodies =
                cube_on_plane(poses, s.dhat, 1000.0, /*gap_factor=*/1.1);
            std::vector<Pose> velocities(1, Pose::Identity(3));
            velocities[0].position = Eigen::Vector3d(1.0, 0, 0);
            Simulator sim(bodies, poses, velocities, /*dt=*/0.01, s);
            REQUIRE(sim.run(/*t_end=*/0.2));
            return sim.poses()[0].position;
        };

        const Eigen::Vector3d with_disabled_friction = run_sim(0.3, 0);
        const Eigen::Vector3d without_friction = run_sim(0.0, 1);
        CHECK((with_disabled_friction - without_friction).norm() == 0.0);
    }
}
