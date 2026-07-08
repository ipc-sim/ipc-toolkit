// Finite-difference and policy tests for the rigid augmented Lagrangian
// (manual derivatives — no autodiff).

#include <catch2/catch_all.hpp>

#include <tests/utils.hpp>

#include <finitediff.hpp>

#include <ipc/dynamics/affine/augmented_lagrangian.hpp>
#include <ipc/dynamics/rigid/augmented_lagrangian.hpp>

using namespace ipc;
using namespace ipc::rigid;

namespace {

std::shared_ptr<RigidBodies> kinematic_cube(std::vector<Pose>& initial_poses)
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    REQUIRE(tests::load_mesh("cube.ply", V, E, F));

    // Non-uniform scaling so J has distinct entries.
    V.col(0).array() *= 0.5;
    V.col(2).array() *= 2.0;

    initial_poses.assign(1, Pose::Identity(3));
    auto bodies = RigidBodies::build_from_meshes(
        { V }, { E }, { F }, { 1000.0 }, initial_poses);
    (*bodies)[0].set_type(RigidBody::Type::KINEMATIC);
    return bodies;
}

} // namespace

TEST_CASE(
    "Rigid augmented Lagrangian gradient and hessian",
    "[rigid][augmented_lagrangian][gradient][hessian]")
{
    std::srand(0);

    std::vector<Pose> initial_poses;
    auto bodies = kinematic_cube(initial_poses);

    // Unreachable satisfied threshold and always-update stall threshold so
    // update() takes the multiplier branch (making λ and Λ nonzero).
    AugmentedLagrangian::Params params;
    params.satisfied_progress = 2.0;
    params.stall_progress = -10.0;
    params.initial_penalty = GENERATE(1e3, 1e6);

    AugmentedLagrangian al(params);

    Pose target = initial_poses[0];
    target.position += Eigen::Vector3d(0.1, -0.2, 0.3);
    target.rotation = Eigen::Vector3d(0.4, -0.1, 0.2);

    const Eigen::VectorXd x0 = Pose::from_poses(initial_poses);
    al.init(*bodies, x0, { target });
    REQUIRE(al.active());

    // Make the multipliers nonzero (multiplier-update branch).
    Eigen::VectorXd x_mid = x0;
    x_mid.head<3>() += Eigen::Vector3d(0.05, -0.1, 0.15);
    x_mid.tail<3>() += Eigen::Vector3d(0.2, -0.05, 0.1);
    al.update(*bodies, x_mid);
    REQUIRE(!al.linear_multipliers()[0].isZero());
    REQUIRE(!al.angular_multipliers()[0].isZero());

    // FD check at a third state.
    Eigen::VectorXd x = x0;
    x.head<3>() += 0.1 * Eigen::Vector3d::Random();
    x.tail<3>() += 0.5 * Eigen::Vector3d::Random();

    const Eigen::VectorXd grad = al.gradient(*bodies, x);
    Eigen::VectorXd fd_grad;
    fd::finite_gradient(
        x, [&](const Eigen::VectorXd& y) { return al(*bodies, y); }, fd_grad);
    CHECK(fd::compare_gradient(grad, fd_grad));

    const Eigen::MatrixXd hess = al.hessian(*bodies, x);
    Eigen::MatrixXd fd_hess;
    fd::finite_jacobian(
        x, [&](const Eigen::VectorXd& y) { return al.gradient(*bodies, y); },
        fd_hess);
    CHECK(fd::compare_hessian(hess, fd_hess));
}

TEST_CASE(
    "Rigid augmented Lagrangian policy", "[rigid][augmented_lagrangian]")
{
    std::vector<Pose> initial_poses;
    auto bodies = kinematic_cube(initial_poses);

    AugmentedLagrangian al((AugmentedLagrangian::Params()));

    Pose target = initial_poses[0];
    target.position += Eigen::Vector3d(1, 0, 0);
    target.rotation = Eigen::Vector3d(0, 0.5, 0);

    const Eigen::VectorXd x0 = Pose::from_poses(initial_poses);
    al.init(*bodies, x0, { target });
    REQUIRE(al.active());
    CHECK(al.linear_penalty() == 1e3);

    // η at the start is 0.
    CHECK(al.linear_progress(*bodies, x0) == Catch::Approx(0.0));
    CHECK(al.angular_progress(*bodies, x0) == Catch::Approx(0.0));

    // η at the target is 1 → both channels satisfy.
    const Eigen::VectorXd x_target =
        Pose::from_poses(std::vector<Pose> { target });
    CHECK(al.linear_progress(*bodies, x_target) == Catch::Approx(1.0));
    CHECK(al.angular_progress(*bodies, x_target) == Catch::Approx(1.0));
    al.update(*bodies, x_target);
    CHECK(al.linear_satisfied());
    CHECK(al.angular_satisfied());
    CHECK(!al.active());

    // Once inactive, the energy is identically zero.
    CHECK(al(*bodies, x0) == 0.0);
    CHECK(al.gradient(*bodies, x0).isZero());

    // Re-init resets; no progress → the penalty doubles up to the cap.
    al.init(*bodies, x0, { target });
    REQUIRE(al.active());
    double prev_kappa = al.linear_penalty();
    for (int i = 0; i < 25; ++i) {
        al.update(*bodies, x0); // η = 0 < stall_progress
        // The doubling stops once κ ≥ max_penalty (checked before doubling,
        // as in the original, so κ can land at up to 2×max).
        CHECK(al.linear_penalty() <= 2e8);
        CHECK(al.linear_penalty() >= prev_kappa);
        prev_kappa = al.linear_penalty();
    }
    const double capped = al.linear_penalty();
    CHECK(capped >= 1e8);
    al.update(*bodies, x0);
    CHECK(al.linear_penalty() == capped); // no further doubling

    // Non-kinematic bodies contribute nothing.
    (*bodies)[0].set_type(RigidBody::Type::DYNAMIC);
    al.init(*bodies, x0, { target });
    CHECK(!al.active());
    CHECK(al(*bodies, x0) == 0.0);
}

TEST_CASE(
    "Affine augmented Lagrangian matches rigid on rotations",
    "[affine][augmented_lagrangian][gradient][hessian]")
{
    std::srand(0);

    std::vector<Pose> initial_poses;
    auto bodies = kinematic_cube(initial_poses);

    // Zero multipliers (no update calls): the affine AL restricted to
    // y = [p; vec(R(θ))] equals the rigid AL since
    // (a − â)ᵀ(J ⊗ I₃)(a − â) = tr((Q − Q̂) J (Q − Q̂)ᵀ).
    rigid::AugmentedLagrangian rigid_al((rigid::AugmentedLagrangian::Params()));
    affine::AugmentedLagrangian affine_al(
        (affine::AugmentedLagrangian::Params()));

    Pose target = initial_poses[0];
    target.position += Eigen::Vector3d(0.3, 0.1, -0.2);
    target.rotation = Eigen::Vector3d(-0.2, 0.4, 0.1);

    const Eigen::VectorXd x0_rigid = Pose::from_poses(initial_poses);
    Eigen::VectorXd x0_affine(12);
    x0_affine.head<3>() = initial_poses[0].position;
    x0_affine.tail<9>() = initial_poses[0].rotation_matrix().reshaped();

    rigid_al.init(*bodies, x0_rigid, { target });
    affine_al.init(*bodies, x0_affine, { target });

    // Evaluate both at a rotation state.
    Pose state = initial_poses[0];
    state.position += Eigen::Vector3d(0.1, -0.05, 0.2);
    state.rotation = Eigen::Vector3d(0.3, -0.2, 0.5);

    const Eigen::VectorXd x_rigid =
        Pose::from_poses(std::vector<Pose> { state });
    Eigen::VectorXd x_affine(12);
    x_affine.head<3>() = state.position;
    x_affine.tail<9>() = state.rotation_matrix().reshaped();

    CHECK(
        rigid_al(*bodies, x_rigid)
        == Catch::Approx(affine_al(*bodies, x_affine)));

    // FD check of the affine AL (fully quadratic).
    const Eigen::VectorXd grad = affine_al.gradient(*bodies, x_affine);
    Eigen::VectorXd fd_grad;
    fd::finite_gradient(
        x_affine,
        [&](const Eigen::VectorXd& y) { return affine_al(*bodies, y); },
        fd_grad);
    CHECK(fd::compare_gradient(grad, fd_grad));

    const Eigen::MatrixXd hess = affine_al.hessian(*bodies, x_affine);
    Eigen::MatrixXd fd_hess;
    fd::finite_jacobian(
        x_affine,
        [&](const Eigen::VectorXd& y) {
            return affine_al.gradient(*bodies, y);
        },
        fd_hess);
    CHECK(fd::compare_hessian(hess, fd_hess));
}
