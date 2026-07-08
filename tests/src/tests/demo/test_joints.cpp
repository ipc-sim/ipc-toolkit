#include <Eigen/Core>
#include <catch2/catch_all.hpp>

#include <tests/utils.hpp>

#include <ipc/dynamics/affine/affine_dof.hpp>
#include <ipc/dynamics/affine/joints.hpp>
#include <simulator.hpp>
#include <ipc/dynamics/rigid/rigid_bodies.hpp>

#include <igl/PI.h>

using namespace ipc;
using namespace ipc::affine;
using BodyDynamics = ipc::demo::Simulator::BodyDynamics;
using ipc::demo::Simulator;

namespace {
/// World position of a world-frame anchor (defined at the initial pose) at
/// the current affine DOFs.
Eigen::Vector3d anchor_position(
    const rigid::RigidBodies& bodies,
    const std::vector<rigid::Pose>& initial_poses,
    Eigen::ConstRef<Eigen::VectorXd> x,
    const size_t body,
    Eigen::ConstRef<Eigen::Vector3d> world_anchor)
{
    const rigid::Pose& p0 = initial_poses[body];
    const Eigen::Vector3d material =
        p0.rotation_matrix().transpose() * (world_anchor - p0.position);
    const Eigen::Vector3d p = x.segment<3>(12 * body);
    const Eigen::Matrix3d A = x.segment<9>(12 * body + 3).reshaped(3, 3);
    return A * material + p;
}

Eigen::VectorXd poses_to_dof(const std::vector<affine::Pose>& poses)
{
    Eigen::VectorXd x(12 * poses.size());
    for (size_t i = 0; i < poses.size(); i++) {
        x.segment<3>(12 * i) = poses[i].position;
        x.segment<9>(12 * i + 3) = poses[i].rotation.reshaped();
    }
    return x;
}
} // namespace

TEST_CASE("Joints require affine dynamics", "[affine][joints][simulator]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    tests::load_mesh("cube.ply", V, E, F);

    std::vector<rigid::Pose> initial_poses(1, rigid::Pose::Identity(3));
    auto bodies = rigid::RigidBodies::build_from_meshes(
        { V }, { E }, { F }, { 1000.0 }, initial_poses);

    auto joints = std::make_shared<JointConstraints>(bodies, initial_poses);
    joints->add_fixed_point(0, Eigen::Vector3d(-0.5, 0, 0));

    // Material-point constraints are nonlinear in the rigid rotation vectors
    Simulator::Settings settings;
    settings.body_dynamics = BodyDynamics::RIGID;
    CHECK_THROWS_AS(
        Simulator(bodies, initial_poses, joints, /*dt=*/0.01, settings),
        std::invalid_argument);
}

TEST_CASE("Joint change of variables", "[affine][joints]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    tests::load_mesh("cube.ply", V, E, F);

    std::vector<rigid::Pose> initial_poses(2, rigid::Pose::Identity(3));
    initial_poses[1].position = Eigen::Vector3d(1, 0, 0);
    auto bodies = rigid::RigidBodies::build_from_meshes(
        { V, V }, { E, E }, { F, F }, { 1000.0, 1000.0 }, initial_poses);

    JointConstraints joints(bodies, initial_poses);
    joints.add_point_connection(0, 1, Eigen::Vector3d(0.5, 0, 0));
    joints.add_sliding_plane(
        1, Eigen::Vector3d(1, 0, 0), Eigen::Vector3d(0, 0, 1));
    REQUIRE(joints.num_constraints() == 4);
    joints.finalize();

    const int n = 24;
    const int m = int(joints.num_constraints());

    // U V = I (round trip)
    for (int trial = 0; trial < 5; trial++) {
        const Eigen::VectorXd x = Eigen::VectorXd::Random(n);
        CHECK(
            (joints.to_full(joints.to_reduced(x)) - x).norm()
            == Catch::Approx(0).margin(1e-12));
    }

    // Any z with pinned head satisfies the constraints exactly
    for (int trial = 0; trial < 5; trial++) {
        Eigen::VectorXd z = Eigen::VectorXd::Random(n);
        z.head(m) = joints.rhs();
        const Eigen::VectorXd x = joints.to_full(z);
        CHECK(joints.residual(x).norm() == Catch::Approx(0).margin(1e-10));
    }

    // The initial configuration satisfies the constraints
    const Eigen::VectorXd x0 = poses_to_dof([&] {
        std::vector<affine::Pose> poses(2);
        for (int i = 0; i < 2; i++) {
            poses[i].position = initial_poses[i].position;
            poses[i].rotation = initial_poses[i].rotation_matrix();
        }
        return poses;
    }());
    CHECK(joints.residual(x0).norm() == Catch::Approx(0).margin(1e-12));
}

TEST_CASE("Affine pendulum on a fixed point", "[affine][joints][simulator]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    tests::load_mesh("cube.ply", V, E, F);

    // Rod hanging horizontally from one end: it should swing down
    Eigen::MatrixXd V_rod = V;
    V_rod.col(0) *= 2.0;

    std::vector<rigid::Pose> initial_poses(1, rigid::Pose::Identity(3));
    auto bodies = rigid::RigidBodies::build_from_meshes(
        { V_rod }, { E }, { F }, { 1000.0 }, initial_poses);

    const Eigen::Vector3d anchor(-1.0, 0, 0); // rod end

    auto joints = std::make_shared<JointConstraints>(bodies, initial_poses);
    joints->add_fixed_point(0, anchor);

    Simulator::Settings settings; // paper-recommended default stiffness
    settings.body_dynamics = BodyDynamics::AFFINE;

    Simulator sim(bodies, initial_poses, joints, /*dt=*/0.01, settings);

    for (int i = 0; i < 25; i++) {
        REQUIRE(sim.step());

        const Eigen::VectorXd x =
            poses_to_dof({ sim.pose_history().back()[0] });

        // The anchor point stays exactly fixed
        const Eigen::Vector3d p =
            anchor_position(*bodies, initial_poses, x, 0, anchor);
        CHECK((p - anchor).norm() == Catch::Approx(0).margin(1e-8));

        // A stays near SO(3)
        const Eigen::Matrix3d A = x.segment<9>(3).reshaped(3, 3);
        CHECK((A.transpose() * A - Eigen::Matrix3d::Identity()).norm() < 0.01);
    }

    // The free end fell (pendulum swings down)
    const Eigen::VectorXd x_final =
        poses_to_dof({ sim.pose_history().back()[0] });
    const Eigen::Vector3d free_end = anchor_position(
        *bodies, initial_poses, x_final, 0, Eigen::Vector3d(1.0, 0, 0));
    CHECK(free_end.y() < -0.1);
    // ...and stays at distance 2 from the anchor (rigidity)
    CHECK((free_end - anchor).norm() == Catch::Approx(2.0).margin(0.02));
}

TEST_CASE("Affine hinge door", "[affine][joints][simulator]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    tests::load_mesh("cube.ply", V, E, F);

    // Body 0: fixed frame; body 1: door attached along a vertical hinge in
    // the gap between the bodies (they must not start in contact — the
    // barrier requires initial separation). Sideways gravity swings the door
    // about the axis.
    std::vector<rigid::Pose> initial_poses(2, rigid::Pose::Identity(3));
    initial_poses[1].position = Eigen::Vector3d(1.2, 0, 0);

    auto bodies = rigid::RigidBodies::build_from_meshes(
        { V, V }, { E, E }, { F, F }, { 1000.0, 1000.0 }, initial_poses);

    const Eigen::Vector3d hinge_p0(0.6, -0.5, 0);
    const Eigen::Vector3d hinge_p1(0.6, 0.5, 0);

    auto joints = std::make_shared<JointConstraints>(bodies, initial_poses);
    joints->add_fixed_body(0);
    joints->add_hinge(0, 1, hinge_p0, hinge_p1);

    Simulator::Settings settings; // paper-recommended default stiffness
    settings.body_dynamics = BodyDynamics::AFFINE;
    settings.gravity = Eigen::Vector3d(0, 0, -10.0); // push the door sideways

    Simulator sim(bodies, initial_poses, joints, /*dt=*/0.01, settings);

    for (int i = 0; i < 25; i++) {
        REQUIRE(sim.step());

        Eigen::VectorXd x(24);
        for (int b = 0; b < 2; b++) {
            const affine::Pose& pose = sim.pose_history().back()[b];
            x.segment<3>(12 * b) = pose.position;
            x.segment<9>(12 * b + 3) = pose.rotation.reshaped();
        }

        // Body 0 stays fixed
        CHECK(
            (x.head<12>() - poses_to_dof({ [&] {
                 affine::Pose pose;
                 pose.position = initial_poses[0].position;
                 pose.rotation = initial_poses[0].rotation_matrix();
                 return pose;
             }() }))
                .norm()
            == Catch::Approx(0).margin(1e-8));

        // Both hinge anchor pairs coincide
        for (const auto& hp : { hinge_p0, hinge_p1 }) {
            const Eigen::Vector3d p_frame =
                anchor_position(*bodies, initial_poses, x, 0, hp);
            const Eigen::Vector3d p_door =
                anchor_position(*bodies, initial_poses, x, 1, hp);
            CHECK((p_frame - p_door).norm() == Catch::Approx(0).margin(1e-8));
        }
    }

    // The door swung about the hinge: its center moved in z
    const affine::Pose& door = sim.pose_history().back()[1];
    CHECK(std::abs(door.position.z()) > 0.05);
}
