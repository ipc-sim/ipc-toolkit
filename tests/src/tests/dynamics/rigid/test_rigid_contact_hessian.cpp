#include <catch2/catch_test_macros.hpp>

#include <finitediff.hpp>

#include <ipc/collisions/normal/normal_collisions.hpp>
#include <ipc/dynamics/rigid/rigid_bodies.hpp>
#include <ipc/dynamics/rigid/rigid_candidates.hpp>
#include <ipc/dynamics/rigid/rigid_contact_hessian.hpp>
#include <ipc/potentials/barrier_potential.hpp>

#include <iostream>

using namespace ipc;
using namespace ipc::rigid;

namespace {

/// A 3×3 grid of isolated points in the y = 0 plane plus one point above it so
/// the inertia is nonsingular (all vertices are codimensional).
Eigen::MatrixXd point_grid()
{
    Eigen::MatrixXd V(10, 3);
    int r = 0;
    for (int i = -1; i <= 1; ++i) {
        for (int j = -1; j <= 1; ++j) {
            V.row(r++) = Eigen::RowVector3d(0.5 * i, 0, 0.5 * j);
        }
    }
    V.row(r) = Eigen::RowVector3d(0, 0.3, 0);
    return V;
}

} // namespace

TEST_CASE(
    "Per-contact rigid Hessian is the exact Hessian without projection",
    "[rigid][hessian]")
{
    constexpr double dhat = 1e-2;

    // Two identical point grids, one hovering just above the other and
    // slightly rotated, so many point-point pairs sit inside dhat.
    std::vector<Pose> poses(2, Pose::identity(3));
    const Eigen::MatrixXd V_grid = point_grid();
    const Eigen::MatrixXi no_elements(0, 2), no_faces(0, 3);
    auto bodies = RigidBodies::build_from_meshes(
        { V_grid, V_grid }, { no_elements, no_elements },
        { no_faces, no_faces }, { 1000.0, 1000.0 }, poses);

    // build_from_meshes moves each body into its principal frame and stores
    // the compensating transform in the pose, so we perturb that pose rather
    // than overwrite it.
    poses[0] = poses[1];
    poses[0].position += Eigen::Vector3d(0.001, 0.004, 0);
    poses[0].rotation += Eigen::Vector3d(2e-3, -1e-3, 3e-3);

    const BarrierPotential barrier(dhat, /*stiffness=*/1.0);
    // Use the rigid body-pair broad phase, as the simulator does: it is what
    // produces the point-point candidates between codimensional vertices.
    RigidCandidates candidates;
    candidates.build(*bodies, poses, /*inflation_radius=*/dhat);
    NormalCollisions collisions;
    collisions.build(candidates, *bodies, bodies->vertices(poses), dhat);
    REQUIRE(collisions.size() > 0);

    // Energy and gradient in rigid DOFs with the collision set held fixed.
    const auto energy = [&](const Eigen::VectorXd& x) {
        return barrier(
            collisions, *bodies, bodies->vertices(Pose::to_poses(x, 3)));
    };
    const auto gradient = [&](const Eigen::VectorXd& x) -> Eigen::VectorXd {
        const std::vector<Pose> p = Pose::to_poses(x, 3);
        return bodies->to_rigid_dof(
            p, barrier.gradient(collisions, *bodies, bodies->vertices(p)));
    };

    const Eigen::VectorXd x = Pose::from_poses(poses);

    Eigen::VectorXd fd_g;
    fd::finite_gradient(x, energy, fd_g);
    CHECK(fd::compare_gradient(gradient(x), fd_g));

    const Eigen::MatrixXd H = per_contact_rigid_hessian(
        *bodies, poses, collisions, barrier, bodies->vertices(poses),
        /*gn_scale=*/1.0, /*geom_scale=*/1.0, PSDProjectionMethod::NONE);
    Eigen::MatrixXd fd_H;
    fd::finite_jacobian(x, gradient, fd_H);

    CHECK(fd::compare_hessian(H, fd_H, 1e-3));
    if (!fd::compare_hessian(H, fd_H, 1e-3)) {
        std::cout << "analytic:\n" << H << "\n\nnumerical:\n" << fd_H << "\n";
    }
}
