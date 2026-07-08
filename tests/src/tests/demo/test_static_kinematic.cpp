// Scenario tests for STATIC and KINEMATIC bodies (DOF masking + augmented
// Lagrangian), run in both RIGID and AFFINE modes.

#include <Eigen/Core>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <tests/utils.hpp>

#include <ipc/dynamics/affine/joints.hpp>
#include <ipc/dynamics/rigid/rigid_bodies.hpp>
#include <ipc/ipc.hpp> // has_intersections
#include <simulator.hpp>

using namespace ipc;
using namespace ipc::rigid;
using ipc::demo::Simulator;

namespace {

Simulator::BodyDynamics generate_mode()
{
    return GENERATE(
        Simulator::BodyDynamics::RIGID, Simulator::BodyDynamics::AFFINE);
}

void load_cube(Eigen::MatrixXd& V, Eigen::MatrixXi& E, Eigen::MatrixXi& F)
{
    REQUIRE(tests::load_mesh("cube.ply", V, E, F));
}

} // namespace

TEST_CASE("Static body support", "[static][simulator]")
{
    const auto mode = generate_mode();
    CAPTURE(int(mode));

    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    load_cube(V, E, F);

    // A static floor cube and a dynamic cube dropped from slightly above.
    std::vector<Pose> initial_poses(2, Pose::Identity(3));
    initial_poses[1].position = Eigen::Vector3d(0.1, 1.1, 0.05);

    auto bodies = RigidBodies::build_from_meshes(
        { V, V }, { E, E }, { F, F }, { 1000.0, 1000.0 }, initial_poses);
    (*bodies)[0].set_type(RigidBody::Type::STATIC);

    Simulator::Settings settings;
    settings.body_dynamics = mode;

    Simulator sim(bodies, initial_poses, /*dt=*/0.01, settings);

    const affine::Pose floor_pose_0 = sim.poses()[0];

    REQUIRE(sim.run(/*t_end=*/0.5));

    // The static floor never moves (bitwise).
    CHECK(sim.poses()[0].position == floor_pose_0.position);
    CHECK(sim.poses()[0].rotation == floor_pose_0.rotation);

    // The dynamic cube rests on top without intersecting.
    CHECK(sim.poses()[1].position.y() > 0.9);
    CHECK(!has_intersections(
        *bodies, bodies->vertices(sim.rigid_pose_history().back())));
}

TEST_CASE("All bodies static", "[static][simulator]")
{
    const auto mode = generate_mode();
    CAPTURE(int(mode));

    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    load_cube(V, E, F);

    std::vector<Pose> initial_poses(1, Pose::Identity(3));
    auto bodies = RigidBodies::build_from_meshes(
        { V }, { E }, { F }, { 1000.0 }, initial_poses);
    (*bodies)[0].set_type(RigidBody::Type::STATIC);

    Simulator::Settings settings;
    settings.body_dynamics = mode;

    Simulator sim(bodies, initial_poses, /*dt=*/0.01, settings);
    const affine::Pose pose_0 = sim.poses()[0];

    REQUIRE(sim.run(/*t_end=*/0.1));
    CHECK(sim.poses()[0].position == pose_0.position);
    CHECK(sim.poses()[0].rotation == pose_0.rotation);
}

TEST_CASE("Kinematic body velocity driven", "[kinematic][simulator]")
{
    const auto mode = generate_mode();
    CAPTURE(int(mode));

    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    load_cube(V, E, F);

    std::vector<Pose> initial_poses(1, Pose::Identity(3));
    auto bodies = RigidBodies::build_from_meshes(
        { V }, { E }, { F }, { 1000.0 }, initial_poses);
    (*bodies)[0].set_type(RigidBody::Type::KINEMATIC);

    std::vector<Pose> initial_velocities(1, Pose::Identity(3));
    initial_velocities[0].position = Eigen::Vector3d(1.0, 0.5, 0);
    initial_velocities[0].rotation = Eigen::Vector3d(0, 0, 0.3);

    Simulator::Settings settings;
    settings.body_dynamics = mode;
    // Gravity must not affect a kinematic body (no inertia/body forces).

    const double dt = 0.01;
    Simulator sim(bodies, initial_poses, initial_velocities, dt, settings);

    const Eigen::Vector3d p0 = initial_poses[0].position;
    const Eigen::Matrix3d Q0 = initial_poses[0].rotation_matrix();
    const Eigen::Vector3d v = initial_velocities[0].position;
    const Eigen::Vector3d omega = initial_velocities[0].rotation;

    const int n = 20;
    for (int i = 0; i < n; ++i) {
        REQUIRE(sim.step());
    }

    // Constant-velocity drive: pₙ = p₀ + n·h·v (within the AL satisfaction
    // tolerance 1 − 0.999 per step) — and no gravity sag.
    const Eigen::Vector3d p_expected = p0 + n * dt * v;
    CHECK(
        (sim.poses()[0].position - p_expected).norm()
        < 2e-3 * n * dt * v.norm());

    // Rotation follows exp(n h [ω]×) Q₀.
    const Eigen::Matrix3d Q_expected =
        Eigen::AngleAxisd(n * dt * omega.norm(), omega.normalized())
            .toRotationMatrix()
        * Q0;
    CHECK((sim.poses()[0].rotation - Q_expected).norm() < 5e-3);
}

TEST_CASE("Kinematic body scripted poses", "[kinematic][simulator]")
{
    const auto mode = generate_mode();
    CAPTURE(int(mode));

    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    load_cube(V, E, F);

    std::vector<Pose> initial_poses(1, Pose::Identity(3));
    auto bodies = RigidBodies::build_from_meshes(
        { V }, { E }, { F }, { 1000.0 }, initial_poses);
    (*bodies)[0].set_type(RigidBody::Type::KINEMATIC);

    // Scripted sinusoidal path (absolute poses, one per step).
    // NOTE: build_from_meshes folds the principal-axes rotation R₀ into
    // initial_poses; compose the script with the post-build pose.
    const int n = 10;
    const double dt = 0.01;
    std::deque<Pose> script;
    for (int i = 1; i <= n; ++i) {
        Pose pose = initial_poses[0];
        pose.position += Eigen::Vector3d(
            0.1 * i * dt, 0.05 * std::sin(0.5 * i), 0);
        script.push_back(pose);
    }
    Simulator::Settings settings;
    settings.body_dynamics = mode;

    Simulator sim(bodies, initial_poses, dt, settings);
    sim.set_kinematic_driver(0, demo::KinematicDriver::scripted(script));

    for (int i = 0; i < n; ++i) {
        const Pose& target = script[i];
        REQUIRE(sim.step());
        // Within the AL satisfaction tolerance of the scripted pose.
        CHECK(
            (sim.poses()[0].position - target.position).norm()
            < 1e-3 * std::max(1.0, target.position.norm()));
    }
}

TEST_CASE("Kinematic body pushes dynamic body", "[kinematic][simulator]")
{
    const auto mode = generate_mode();
    CAPTURE(int(mode));

    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    load_cube(V, E, F);

    // Kinematic cube slides toward a dynamic cube resting to its right.
    std::vector<Pose> initial_poses(2, Pose::Identity(3));
    initial_poses[1].position = Eigen::Vector3d(1.2, 0, 0);

    auto bodies = RigidBodies::build_from_meshes(
        { V, V }, { E, E }, { F, F }, { 1000.0, 1000.0 }, initial_poses);
    (*bodies)[0].set_type(RigidBody::Type::KINEMATIC);

    std::vector<Pose> initial_velocities(2, Pose::Identity(3));
    initial_velocities[0].position = Eigen::Vector3d(1.0, 0, 0);

    Simulator::Settings settings;
    settings.body_dynamics = mode;
    settings.gravity.setZero(); // isolate the push

    const double dt = 0.01;
    Simulator sim(bodies, initial_poses, initial_velocities, dt, settings);

    const double x1_start = initial_poses[1].position.x();
    for (int i = 0; i < 50; ++i) {
        CAPTURE(i);
        REQUIRE(sim.step());
        CHECK(!has_intersections(
            *bodies, bodies->vertices(sim.rigid_pose_history().back())));
    }

    // The kinematic cube advanced ~0.5 m; the gap was 0.2 m, so the dynamic
    // cube must have been pushed forward without tunneling.
    CHECK(sim.poses()[0].position.x() > 0.45);
    CHECK(sim.poses()[1].position.x() > x1_start + 0.2);
    CHECK(
        sim.poses()[1].position.x() - sim.poses()[0].position.x()
        > 0.99); // still at least one cube width apart (no tunneling)
}

TEST_CASE("Partial DOF fixing", "[static][simulator]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    load_cube(V, E, F);

    // Fix the rotation of a cube; push it off-center so an unconstrained
    // cube would spin.
    std::vector<Pose> initial_poses(1, Pose::Identity(3));
    VectorMax6b is_dof_fixed = VectorMax6b::Zero(6);
    is_dof_fixed.tail<3>().setOnes(); // fix all rotation

    auto bodies = RigidBodies::build_from_meshes(
        { V }, { E }, { F }, { 1000.0 }, initial_poses,
        /*convert_planes=*/false, { is_dof_fixed });

    // Off-center force: a torque plus a linear force.
    (*bodies)[0].set_external_force(
        Pose(Eigen::Vector3d(500.0, 0, 0), Eigen::Vector3d(0, 0, 800.0)));

    const auto mode = generate_mode();
    CAPTURE(int(mode));
    Simulator::Settings settings;
    settings.body_dynamics = mode;
    settings.gravity.setZero();

    Simulator sim(bodies, initial_poses, /*dt=*/0.01, settings);

    const affine::Pose pose_0 = sim.poses()[0];
    REQUIRE(sim.run(/*t_end=*/0.2));

    // Orientation unchanged; position responded to the force.
    CHECK(
        (sim.poses()[0].rotation - pose_0.rotation).norm()
        == Catch::Approx(0).margin(1e-9));
    CHECK(sim.poses()[0].position.x() > pose_0.position.x() + 1e-4);
}

TEST_CASE("Kinematic max time expiry", "[kinematic][simulator]")
{
    const auto mode = generate_mode();
    CAPTURE(int(mode));

    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    load_cube(V, E, F);

    std::vector<Pose> initial_poses(1, Pose::Identity(3));
    auto bodies = RigidBodies::build_from_meshes(
        { V }, { E }, { F }, { 1000.0 }, initial_poses);
    (*bodies)[0].set_type(RigidBody::Type::KINEMATIC);

    const double dt = 0.01;

    std::vector<Pose> initial_velocities(1, Pose::Identity(3));
    initial_velocities[0].position = Eigen::Vector3d(1.0, 0, 0);

    Simulator::Settings settings;
    settings.body_dynamics = mode;
    settings.gravity.setZero();

    Simulator sim(bodies, initial_poses, initial_velocities, dt, settings);
    sim.set_kinematic_driver(
        0, demo::KinematicDriver::velocity_driven(/*max_time=*/2 * dt));

    REQUIRE(sim.step()); // kinematic step 1
    REQUIRE(sim.step()); // kinematic step 2 (max time reached)
    CHECK((*bodies)[0].type() == RigidBody::Type::KINEMATIC);

    REQUIRE(sim.step()); // converts to STATIC at step start
    CHECK((*bodies)[0].type() == RigidBody::Type::STATIC);

    const affine::Pose pose_frozen = sim.poses()[0];
    REQUIRE(sim.step());
    CHECK(sim.poses()[0].position == pose_frozen.position);
    CHECK(sim.poses()[0].rotation == pose_frozen.rotation);

    // reset() restores the KINEMATIC type (and re-arms the driver, so the
    // body drives again and expires a second time).
    sim.reset();
    CHECK((*bodies)[0].type() == RigidBody::Type::KINEMATIC);
    REQUIRE(sim.step());
    REQUIRE(sim.step());
    CHECK((*bodies)[0].type() == RigidBody::Type::KINEMATIC);
    REQUIRE(sim.step());
    CHECK((*bodies)[0].type() == RigidBody::Type::STATIC);
}

TEST_CASE("Joints reject non-dynamic bodies", "[kinematic][joints][simulator]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    load_cube(V, E, F);

    std::vector<Pose> initial_poses(2, Pose::Identity(3));
    initial_poses[1].position = Eigen::Vector3d(1.5, 0, 0);

    auto bodies = RigidBodies::build_from_meshes(
        { V, V }, { E, E }, { F, F }, { 1000.0, 1000.0 }, initial_poses);
    (*bodies)[0].set_type(RigidBody::Type::KINEMATIC);

    auto joints =
        std::make_shared<affine::JointConstraints>(bodies, initial_poses);
    joints->add_point_connection(0, 1, Eigen::Vector3d(0.75, 0, 0));

    Simulator::Settings settings;
    settings.body_dynamics = Simulator::BodyDynamics::AFFINE;

    CHECK_THROWS_AS(
        Simulator(bodies, initial_poses, joints, /*dt=*/0.01, settings),
        std::invalid_argument);
}
