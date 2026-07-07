#include <Eigen/Core>
#include <catch2/catch_all.hpp>

#include <tests/utils.hpp>

#include <ipc/dynamics/affine/affine_dof.hpp>
#include <ipc/dynamics/affine/inertial_term.hpp>
#include <ipc/dynamics/time_integration/bdf.hpp>
#include <simulator.hpp>
#include <ipc/dynamics/rigid/inertial_term.hpp>
#include <ipc/dynamics/rigid/rigid_bodies.hpp>
#include <ipc/potentials/barrier_potential.hpp>

#include <finitediff.hpp>
#include <igl/PI.h>

#include <iostream>

using namespace ipc;
using namespace ipc::affine;
using BodyDynamics = ipc::demo::Simulator::BodyDynamics;
using ipc::demo::Simulator;

namespace {
/// Embed rigid poses (6 DOF) into affine DOFs (12 DOF, column-major vec(A)).
Eigen::VectorXd rigid_to_affine_dof(Eigen::ConstRef<Eigen::VectorXd> x_rigid)
{
    const size_t num_bodies = x_rigid.size() / 6;
    Eigen::VectorXd x(12 * num_bodies);
    for (size_t i = 0; i < num_bodies; ++i) {
        x.segment<3>(12 * i) = x_rigid.segment<3>(6 * i);
        x.segment<9>(12 * i + 3) =
            rigid::rotation_vector_to_matrix(x_rigid.segment<3>(6 * i + 3))
                .reshaped();
    }
    return x;
}
} // namespace

TEST_CASE(
    "Affine inertial term matches rigid on rigid states", "[affine][inertial]")
{
    // FD tests cannot catch a wrong-but-symmetric mass matrix; instead check
    // the analytic identity: for x embedding a rigid state (A ∈ SO(3)), the
    // affine inertial energy ½(x−x̂)ᵀM(x−x̂) must equal the rigid inertial
    // energy ½m‖p−p̂‖² + ½tr((Q−Q̂)J(Q−Q̂)ᵀ) when both share the same
    // predicted state.
    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    tests::load_mesh("cube.ply", V, E, F);

    // Stretch the cube so the moment of inertia is anisotropic (this catches
    // convention errors that a cube's J ∝ I would hide)
    Eigen::MatrixXd V_box = V;
    V_box.col(0) *= 2.0;
    V_box.col(1) *= 0.5;

    std::vector<rigid::Pose> initial_poses(2, rigid::Pose::Identity(3));
    initial_poses[1].position = Eigen::Vector3d(3, 0, 0);
    auto bodies = rigid::RigidBodies::build_from_meshes(
        { V_box, V }, { E, E }, { F, F }, { 1000.0, 500.0 }, initial_poses);

    // Random rigid state and velocity for the shared integrator
    const Eigen::VectorXd x0_rigid = Eigen::VectorXd::Random(12);
    const Eigen::VectorXd x0 = rigid_to_affine_dof(x0_rigid);
    const Eigen::VectorXd v0 = 0.1 * Eigen::VectorXd::Random(24);
    const Eigen::VectorXd a0 = Eigen::VectorXd::Zero(24);

    auto integrator = std::make_shared<dynamics::BDF>(/*order=*/1);
    integrator->set_dt(0.01);
    integrator->init(x0, v0, a0, /*num_bodies=*/2);

    affine::InertialTerm affine_inertial(*bodies, integrator);
    affine_inertial.update(*bodies);

    rigid::InertialTerm rigid_inertial(integrator);
    rigid_inertial.update(*bodies);

    for (int trial = 0; trial < 10; ++trial) {
        const Eigen::VectorXd x_rigid = Eigen::VectorXd::Random(12);
        const Eigen::VectorXd x = rigid_to_affine_dof(x_rigid);

        const double affine_energy = affine_inertial(*bodies, x);
        const double rigid_energy = rigid_inertial(*bodies, x_rigid);

        CHECK(
            affine_energy
            == Catch::Approx(rigid_energy).epsilon(1e-10).margin(1e-12));
    }
}

TEST_CASE("Affine DOF jacobian", "[affine][dof]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    tests::load_mesh("cube.ply", V, E, F);

    std::vector<rigid::Pose> initial_poses(2, rigid::Pose::Identity(3));
    initial_poses[1].position = Eigen::Vector3d(2, 0, 0);
    auto bodies = rigid::RigidBodies::build_from_meshes(
        { V, V }, { E, E }, { F, F }, { 1000.0, 1000.0 }, initial_poses);

    const Eigen::SparseMatrix<double> J = affine_jacobian(*bodies);

    const Eigen::VectorXd x = Eigen::VectorXd::Random(24);

    // V(x) is linear, so the FD jacobian must match exactly
    Eigen::MatrixXd J_fd;
    fd::finite_jacobian(
        x,
        [&](const Eigen::VectorXd& x_) -> Eigen::VectorXd {
            return vertices(*bodies, x_).reshaped<Eigen::RowMajor>();
        },
        J_fd);

    CHECK(fd::compare_jacobian(Eigen::MatrixXd(J), J_fd));
}

TEST_CASE("Affine simulator gradient and hessian", "[affine][simulator]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    tests::load_mesh("cube.ply", V, E, F);

    // Two cubes close enough that the barrier is active
    std::vector<rigid::Pose> initial_poses(2, rigid::Pose::Identity(3));
    initial_poses[1].position = Eigen::Vector3d(1.0005, 0, 0);

    auto bodies = rigid::RigidBodies::build_from_meshes(
        { V, V }, { E, E }, { F, F }, { 1000.0, 1000.0 }, initial_poses);

    // Use a mild orthogonality stiffness: derivative correctness is
    // scale-independent, and the paper-recommended κ dwarfs the other terms,
    // hurting the finite-difference conditioning.
    Simulator::Settings settings;
    settings.body_dynamics = BodyDynamics::AFFINE;
    settings.orthogonality_stiffness = 1e6;

    Simulator sim(bodies, initial_poses, /*dt=*/0.01, settings);

    // Perturbed affine state (slightly non-rigid A)
    Eigen::VectorXd x(24);
    for (size_t i = 0; i < 2; ++i) {
        x.segment<3>(12 * i) = initial_poses[i].position;
        x.segment<9>(12 * i + 3) =
            initial_poses[i].rotation_matrix().reshaped();
    }
    // Fixed seed: an unlucky perturbation can separate the cubes beyond dhat,
    // deactivating the barrier this test requires.
    std::srand(0);
    x += 0.001 * Eigen::VectorXd::Random(24);

    sim.initialize_step();
    sim.update_collisions(x);
    REQUIRE(sim.normal_collisions().size() > 0);

    // --- Gradient -----------------------------------------------------------

    Eigen::VectorXd g;
    sim.gradient(x, g);

    Eigen::VectorXd fd_g;
    fd::finite_gradient(
        x, [&](const Eigen::VectorXd& x_) { return sim.value(x_); }, fd_g);

    CHECK(fd::compare_gradient(g, fd_g));
    if (!fd::compare_gradient(g, fd_g)) {
        std::cout << "analytic:\n" << g.transpose() << "\n\n";
        std::cout << "numerical:\n" << fd_g.transpose() << "\n\n";
    }

    // --- Hessian
    // --------------------------------------------------------------

    sim.set_project_to_psd(false); // FD comparison needs the exact Hessian
    Eigen::SparseMatrix<double> H_sparse;
    sim.hessian(x, H_sparse);
    const Eigen::MatrixXd H = H_sparse;

    Eigen::MatrixXd fd_H;
    fd::finite_jacobian(
        x,
        [&](const Eigen::VectorXd& x_) {
            Eigen::VectorXd g_;
            sim.gradient(x_, g_);
            return g_;
        },
        fd_H);

    CHECK(fd::compare_hessian(H, fd_H));
    if (!fd::compare_hessian(H, fd_H)) {
        std::cout << "analytic:\n" << H << "\n\n";
        std::cout << "numerical:\n" << fd_H << "\n\n";
    }
}

TEST_CASE("Affine simulator cube drop", "[affine][simulator]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    tests::load_mesh("cube.ply", V, E, F);

    std::vector<rigid::Pose> initial_poses(1, rigid::Pose::Identity(3));
    initial_poses[0].position = Eigen::Vector3d(0, 0.6, 0);

    auto bodies = rigid::RigidBodies::build_from_meshes(
        { V }, { E }, { F }, { 1000.0 }, initial_poses);
    bodies->planes.emplace_back(Eigen::Vector3d(0, 1, 0), 0); // ground y = 0

    Simulator::Settings settings; // paper-recommended default stiffness
    settings.body_dynamics = BodyDynamics::AFFINE;

    Simulator sim(bodies, initial_poses, /*dt=*/0.01, settings);

    for (int i = 0; i < 25; ++i) {
        REQUIRE(sim.step());

        const rigid::AffinePose& pose = sim.pose_history().back()[0];

        // A stays near SO(3) under the stiff orthogonality potential
        const Eigen::Matrix3d AtA_minus_I =
            pose.rotation.transpose() * pose.rotation
            - Eigen::Matrix3d::Identity();
        CHECK(AtA_minus_I.norm() < 0.01);
    }

    // The cube must rest on (not penetrate) the ground
    const Eigen::MatrixXd V_final = vertices(*bodies, [&] {
        Eigen::VectorXd x(12);
        const rigid::AffinePose& pose = sim.pose_history().back()[0];
        x.head<3>() = pose.position;
        x.tail<9>() = pose.rotation.reshaped();
        return x;
    }());
    CHECK((V_final.col(1).array() > -1e-10).all());
    // It fell and is near the ground
    CHECK(V_final.col(1).minCoeff() < 0.1);
}
