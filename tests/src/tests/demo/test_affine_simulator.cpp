#include <Eigen/Core>
#include <catch2/catch_all.hpp>

#include <tests/utils.hpp>

#include <ipc/dynamics/affine/affine_dof.hpp>
#include <ipc/dynamics/time_integration/bdf.hpp>
#include <simulator.hpp>
#include <ipc/dynamics/rigid/rigid_bodies.hpp>
#include <ipc/potentials/barrier_potential.hpp>

#include <finitediff.hpp>
#include <igl/PI.h>

#include <cmath>
#include <iostream>

using namespace ipc;
using namespace ipc::affine;
using BodyDynamics = ipc::demo::Simulator::BodyDynamics;
using ipc::demo::Simulator;

// NOTE: The inertial-term convention (the affine mass matrix equals the rigid
// trace form on rigid states) is covered by the change-of-variables and body-
// potentials unit tests ([to_affine], [body_potentials]).

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
    // Deterministic pseudo-random perturbation: an unlucky perturbation can
    // separate the cubes beyond dhat, deactivating the barrier this test
    // requires. std::rand() is not used here because it is implementation-
    // defined and seeding it with std::srand(0) does not reproduce the same
    // sequence across compilers/platforms (e.g., MSVC vs. glibc/libc++),
    // which previously made this test flaky on Windows CI.
    Eigen::VectorXd perturbation(24);
    for (int i = 0; i < perturbation.size(); ++i) {
        perturbation(i) = std::sin(12.9898 * (i + 1));
    }
    x += 0.001 * perturbation;

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

        const affine::Pose& pose = sim.pose_history().back()[0];

        // A stays near SO(3) under the stiff orthogonality potential
        const Eigen::Matrix3d AtA_minus_I =
            pose.rotation.transpose() * pose.rotation
            - Eigen::Matrix3d::Identity();
        CHECK(AtA_minus_I.norm() < 0.01);
    }

    // The cube must rest on (not penetrate) the ground
    const Eigen::MatrixXd V_final = vertices(*bodies, [&] {
        Eigen::VectorXd x(12);
        const affine::Pose& pose = sim.pose_history().back()[0];
        x.head<3>() = pose.position;
        x.tail<9>() = pose.rotation.reshaped();
        return x;
    }());
    CHECK((V_final.col(1).array() > -1e-10).all());
    // It fell and is near the ground
    CHECK(V_final.col(1).minCoeff() < 0.1);
}
