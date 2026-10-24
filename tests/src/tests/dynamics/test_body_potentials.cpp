#include <catch2/catch_all.hpp>

#include <tests/utils.hpp>

#include <finitediff.hpp>

#include <ipc/dynamics/body_potentials.hpp>
#include <ipc/dynamics/to_affine.hpp>
#include <ipc/dynamics/rigid/rigid_bodies.hpp>
#include <ipc/dynamics/time_integration/bdf.hpp>

#include <iostream>
#include <memory>

using namespace ipc;
using namespace ipc::dynamics;

namespace {
/// Embed a rigid state (6 DOF) into affine DOFs (12 DOF, column-major vec(A)).
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
    "BodyPotentials gradient and hessian",
    "[affine][body_potentials][gradient][hessian]")
{
    const bool rigid_cov = GENERATE(true, false);
    CAPTURE(rigid_cov);

    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    REQUIRE(tests::load_mesh("cube.ply", V, E, F));
    // Anisotropic so the inertia tensor has distinct entries.
    V.col(0) *= 2.0;
    V.col(1) *= 0.5;

    std::vector<rigid::Pose> initial_poses(2, rigid::Pose::Identity(3));
    initial_poses[1].position = Eigen::Vector3d(3, 0, 0);
    auto bodies = rigid::RigidBodies::build_from_meshes(
        { V, V }, { E, E }, { F, F }, { 1000.0, 500.0 }, initial_poses);

    // Shared integrator with an affine-shaped state and nonzero velocity.
    const Eigen::VectorXd x0 = rigid_to_affine_dof(Eigen::VectorXd::Random(12));
    const Eigen::VectorXd v0 = 0.1 * Eigen::VectorXd::Random(24);
    const Eigen::VectorXd a0 = Eigen::VectorXd::Zero(24);
    auto integrator = std::make_shared<dynamics::BDF>(/*order=*/1);
    integrator->set_dt(0.01);
    integrator->init(x0, v0, a0, /*num_bodies=*/2);

    std::shared_ptr<dynamics::ToAffine> to_affine;
    if (rigid_cov) {
        to_affine = std::make_shared<dynamics::RigidToAffine>(3, 2);
    } else {
        to_affine = std::make_shared<dynamics::AffineToAffine>(3, 2);
    }

    BodyPotentials potential(
        *bodies, integrator, to_affine, /*orthogonality_stiffness=*/1e3);
    potential.set_gravity(Eigen::Vector3d(0, -9.81, 0));
    potential.update(*bodies);

    const double s = 1e-2; // arbitrary positive scaling

    // Evaluate at a perturbed state (reduced DOFs for rigid, affine for ABD).
    Eigen::VectorXd x;
    if (rigid_cov) {
        x = 0.3 * Eigen::VectorXd::Random(12); // [p; θ] per body
    } else {
        x = rigid_to_affine_dof(0.3 * Eigen::VectorXd::Random(12));
    }

    const auto energy = [&](const Eigen::VectorXd& x_) {
        return potential.energy(*bodies, x_, s);
    };
    const auto gradient = [&](const Eigen::VectorXd& x_) {
        return potential.gradient(*bodies, x_, s);
    };

    // --- Gradient -----------------------------------------------------------
    const Eigen::VectorXd g = gradient(x);
    Eigen::VectorXd fd_g;
    fd::finite_gradient(x, energy, fd_g);
    CHECK(fd::compare_gradient(g, fd_g));
    if (!fd::compare_gradient(g, fd_g)) {
        std::cout << "analytic:\n" << g.transpose() << "\n";
        std::cout << "numerical:\n" << fd_g.transpose() << "\n";
    }

    // --- Hessian (exact, no PSD projection) ---------------------------------
    const Eigen::SparseMatrix<double> H =
        potential.hessian(*bodies, x, s, PSDProjectionMethod::NONE);
    Eigen::MatrixXd fd_H;
    fd::finite_jacobian(x, gradient, fd_H);
    CHECK(fd::compare_hessian(Eigen::MatrixXd(H), fd_H, 1e-3));
    if (!fd::compare_hessian(Eigen::MatrixXd(H), fd_H, 1e-3)) {
        std::cout << "analytic:\n" << Eigen::MatrixXd(H) << "\n\n";
        std::cout << "numerical:\n" << fd_H << "\n";
    }
}
