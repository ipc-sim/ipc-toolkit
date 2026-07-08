// End-to-end 2D simulation tests (rigid and affine).

#include <Eigen/Core>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <tests/utils.hpp>

#include <ipc/collisions/normal/normal_collisions.hpp>
#include <ipc/dynamics/affine/joints.hpp>
#include <ipc/dynamics/rigid/rigid_bodies.hpp>
#include <simulator.hpp>

#include <finitediff.hpp>

#include <iostream>

using namespace ipc;
using namespace ipc::rigid;
using ipc::demo::Simulator;

namespace {

/// A unit square outline (4 vertices, 4 edges, no faces).
void unit_square(Eigen::MatrixXd& V, Eigen::MatrixXi& E, Eigen::MatrixXi& F)
{
    V.resize(4, 2);
    V << -0.5, -0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5;
    E.resize(4, 2);
    E << 0, 1, 1, 2, 2, 3, 3, 0;
    F.resize(0, 3);
}

Simulator::BodyDynamics generate_mode()
{
    return GENERATE(
        Simulator::BodyDynamics::RIGID, Simulator::BodyDynamics::AFFINE);
}

} // namespace

TEST_CASE("2D ballistic motion", "[2d][rigid][simulator]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    unit_square(V, E, F);

    std::vector<Pose> initial_poses(1, Pose::Identity(2));
    initial_poses[0].position = Eigen::Vector2d(0, 100.0); // free flight

    auto bodies = RigidBodies::build_from_meshes(
        { V }, { E }, { F }, { 1000.0 }, initial_poses);

    std::vector<Pose> initial_velocities(1, Pose::Identity(2));
    initial_velocities[0].position = Eigen::Vector2d(1.0, 2.0);
    initial_velocities[0].rotation = VectorMax3d::Constant(1, 0.7); // ω

    Simulator::Settings settings;
    settings.solver_params["grad_norm_tol"] = 1e-10;

    const double dt = 0.01;
    Simulator sim(bodies, initial_poses, initial_velocities, dt, settings);

    const Eigen::Vector2d p0 = initial_poses[0].position;
    const double theta0 = initial_poses[0].rotation(0);
    const Eigen::Vector2d v0 = initial_velocities[0].position;
    const double omega = initial_velocities[0].rotation(0);
    const Eigen::Vector2d g = settings.gravity.head<2>();

    const int n = 10;
    for (int i = 0; i < n; ++i) {
        REQUIRE(sim.step());
    }

    // Implicit Euler: pₙ = p₀ + n dt v₀ + dt² g n(n+1)/2.
    const Eigen::Vector2d pn_expected =
        p0 + n * dt * v0 + dt * dt * g * (n * (n + 1) / 2.0);
    CHECK(
        (Eigen::Vector2d(sim.poses()[0].position) - pn_expected).norm()
        == Catch::Approx(0).margin(1e-6));

    const Pose pose_n = sim.rigid_poses()[0];
    // θₙ ≈ θ₀ + n dt ω. The rotation is integrated in matrix space (as in 3D),
    // so the recovered angle differs from the linear-in-angle value by an
    // O((dt ω)³) per-step amount; check with a correspondingly looser margin.
    const double theta_expected =
        std::remainder(theta0 + n * dt * omega, 2 * M_PI);
    CHECK(pose_n.rotation(0) == Catch::Approx(theta_expected).margin(1e-4));
}

TEST_CASE("2D simulator gradient and hessian", "[2d][simulator]")
{
    const auto mode = generate_mode();
    CAPTURE(int(mode));

    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    unit_square(V, E, F);

    Simulator::Settings settings;
    settings.body_dynamics = mode;

    // Two dynamic squares in contact (0.5·dhat gap). NOTE: both dynamic —
    // a STATIC body's pinned DOFs are zeroed in gradient() by design, which
    // would (correctly) disagree with the FD of value().
    std::vector<Pose> initial_poses(2, Pose::Identity(2));
    initial_poses[1].position = Eigen::Vector2d(0.1, 1.0 + 0.5 * settings.dhat);

    auto bodies = RigidBodies::build_from_meshes(
        { V, V }, { E, E }, { F, F }, { 1000.0, 1000.0 }, initial_poses);

    Simulator sim(bodies, initial_poses, /*dt=*/0.01, settings);

    const Eigen::VectorXd x0 = sim.kinematics().dof(initial_poses);
    sim.initialize_step();
    sim.update_collisions(x0);
    REQUIRE(!sim.normal_collisions().empty());

    // Perturb the dynamic body (translation + rotation) to exercise the 2D
    // chain rule (including the rigid second-order term).
    Eigen::VectorXd x = x0;
    const int body_ndof = int(x.size()) / 2;
    x(body_ndof) += 1e-4;
    x(body_ndof + int(bodies->dim())) += 1e-3;

    Eigen::VectorXd grad;
    sim.gradient(x, grad);

    Eigen::VectorXd fd_grad;
    fd::finite_gradient(
        x, [&](const Eigen::VectorXd& y) { return sim.value(y); }, fd_grad);
    CHECK(fd::compare_gradient(grad, fd_grad));
    if (!fd::compare_gradient(grad, fd_grad)) {
        std::cout << "analytic:\n" << grad << "\n\nnumeric:\n" << fd_grad;
    }

    sim.set_project_to_psd(false);
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
    CHECK(fd::compare_hessian(hess, fd_hess, 1e-2));
}

TEST_CASE("2D box drop", "[2d][simulator]")
{
    const auto mode = generate_mode();
    CAPTURE(int(mode));

    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    unit_square(V, E, F);

    // A static floor square and a dynamic square dropped from above.
    std::vector<Pose> initial_poses(2, Pose::Identity(2));
    initial_poses[1].position = Eigen::Vector2d(0.1, 1.3);

    auto bodies = RigidBodies::build_from_meshes(
        { V, V }, { E, E }, { F, F }, { 1000.0, 1000.0 }, initial_poses);
    (*bodies)[0].set_type(RigidBody::Type::STATIC);

    Simulator::Settings settings;
    settings.body_dynamics = mode;

    Simulator sim(bodies, initial_poses, /*dt=*/0.01, settings);
    REQUIRE(sim.run(/*t_end=*/0.75));

    // The dynamic square rests on top of the floor (no tunneling, no
    // penetration): its center stays about one edge length above the floor.
    CHECK(sim.poses()[1].position(1) > 0.99);
    CHECK(sim.poses()[1].position(1) < 1.1);

    if (mode == Simulator::BodyDynamics::AFFINE) {
        // A stays (numerically) a rotation.
        const MatrixMax3d A = sim.poses()[1].rotation;
        CHECK(
            (A.transpose() * A - Eigen::Matrix2d::Identity()).norm()
            == Catch::Approx(0).margin(1e-4));
    }
}

TEST_CASE("2D friction slide", "[2d][friction][simulator]")
{
    const double mu = 0.3, g = 9.81, v0 = 1.0;

    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    unit_square(V, E, F);

    Simulator::Settings settings;
    settings.friction_coefficient = mu;
    settings.friction_iterations = -1;
    settings.velocity_conv_tol = 1e-3;

    // Wide static floor (three squares would be nicer, but a stretched one
    // works): scale the floor square horizontally.
    Eigen::MatrixXd V_floor = V;
    V_floor.col(0) *= 10.0;

    std::vector<Pose> initial_poses(2, Pose::Identity(2));
    initial_poses[1].position = Eigen::Vector2d(-2, 1.0 + 1.1 * settings.dhat);

    auto bodies = RigidBodies::build_from_meshes(
        { V_floor, V }, { E, E }, { F, F }, { 1000.0, 1000.0 }, initial_poses);
    (*bodies)[0].set_type(RigidBody::Type::STATIC);

    std::vector<Pose> initial_velocities(2, Pose::Identity(2));
    initial_velocities[1].position = Eigen::Vector2d(v0, 0);

    const double dt = 0.01;
    Simulator sim(bodies, initial_poses, initial_velocities, dt, settings);

    const double x_start = initial_poses[1].position(0);
    const double stopping_time = v0 / (mu * g);
    const double stopping_dist = v0 * v0 / (2 * mu * g);

    REQUIRE(sim.run(/*t_end=*/2 * stopping_time));
    const double slide = sim.poses()[1].position(0) - x_start;
    CHECK(slide == Catch::Approx(stopping_dist).margin(3 * v0 * dt));
}

TEST_CASE("2D kinematic body", "[2d][kinematic][simulator]")
{
    const auto mode = generate_mode();
    CAPTURE(int(mode));

    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    unit_square(V, E, F);

    std::vector<Pose> initial_poses(1, Pose::Identity(2));
    auto bodies = RigidBodies::build_from_meshes(
        { V }, { E }, { F }, { 1000.0 }, initial_poses);
    (*bodies)[0].set_type(RigidBody::Type::KINEMATIC);

    std::vector<Pose> initial_velocities(1, Pose::Identity(2));
    initial_velocities[0].position = Eigen::Vector2d(1.0, 0.5);
    initial_velocities[0].rotation = VectorMax3d::Constant(1, 0.3);

    Simulator::Settings settings;
    settings.body_dynamics = mode;

    const double dt = 0.01;
    Simulator sim(bodies, initial_poses, initial_velocities, dt, settings);

    const Eigen::Vector2d p0 = initial_poses[0].position;
    const double theta0 = initial_poses[0].rotation(0);

    const int n = 10;
    for (int i = 0; i < n; ++i) {
        REQUIRE(sim.step());
    }

    const Eigen::Vector2d p_expected = p0 + n * dt * Eigen::Vector2d(1.0, 0.5);
    CHECK(
        (Eigen::Vector2d(sim.poses()[0].position) - p_expected).norm()
        < 2e-3 * n * dt);
    CHECK(
        sim.rigid_poses()[0].rotation(0)
        == Catch::Approx(theta0 + n * dt * 0.3).margin(5e-3));
}

TEST_CASE("2D joints", "[2d][affine][joints][simulator]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    unit_square(V, E, F);

    std::vector<Pose> initial_poses(1, Pose::Identity(2));
    auto bodies = RigidBodies::build_from_meshes(
        { V }, { E }, { F }, { 1000.0 }, initial_poses);

    auto joints =
        std::make_shared<affine::JointConstraints>(bodies, initial_poses);

    // A 2D hinge is a point connection.
    CHECK_THROWS_AS(
        joints->add_hinge(
            0, 0, Eigen::Vector2d(0, 0.5), Eigen::Vector2d(0, -0.5)),
        std::invalid_argument);

    // Pendulum: pin the top-left corner of the square; it swings under
    // gravity, keeping the pinned point fixed.
    const Eigen::Vector2d pivot =
        initial_poses[0].rotation_matrix() * Eigen::Vector2d(-0.5, 0.5)
        + Eigen::Vector2d(initial_poses[0].position);
    joints->add_fixed_point(0, pivot);

    Simulator::Settings settings;
    settings.body_dynamics = Simulator::BodyDynamics::AFFINE;

    Simulator sim(bodies, initial_poses, joints, /*dt=*/0.01, settings);

    double max_swing = 0;
    for (int i = 0; i < 50; ++i) {
        REQUIRE(sim.step());
        // The pinned material point stays at the pivot.
        const auto& pose = sim.poses()[0];
        const Eigen::Vector2d p = Eigen::Matrix2d(pose.rotation)
                * (initial_poses[0].rotation_matrix().transpose()
                   * (pivot - Eigen::Vector2d(initial_poses[0].position)))
            + Eigen::Vector2d(pose.position);
        CHECK((p - pivot).norm() == Catch::Approx(0).margin(1e-8));
        max_swing = std::max(
            max_swing,
            std::abs(
                Eigen::Vector2d(pose.position).x()
                - initial_poses[0].position(0)));
    }
    CHECK(max_swing > 0.01); // it actually swings
}
