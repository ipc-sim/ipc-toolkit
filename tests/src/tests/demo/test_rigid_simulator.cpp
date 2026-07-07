// Test the mass utilities.

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <tests/utils.hpp>

#include <ipc/config.hpp>
#include <ipc/ipc.hpp> // has_intersections
#include <simulator.hpp>
#include <ipc/dynamics/rigid/rigid_bodies.hpp>
#ifdef IPC_TOOLKIT_WITH_GLTF
#include <ipc/io/write_gltf.hpp>
#endif
#include <ipc/potentials/barrier_potential.hpp>

#include <finitediff.hpp>
#include <igl/PI.h>

#include <iostream>

using namespace ipc;
using namespace ipc::rigid;
using ipc::demo::Simulator;

namespace ipc::tests {
void simulator_test(Simulator& sim)
{
    Eigen::VectorXd x(18);
    x << -0.31691, 0.358068, 0.0934453, 1.15926, 2.13156, 2.6372, //
        0.32609, 0.41217, 0.256955, 0.515056, 1.56975, 0.399389,  //
        0.00250424, 0.11456, -0.00196945, 0.289872, -2.08816, 2.06641;

    // --- Check barrier potential gradient ------------------------------------

    Eigen::MatrixXd V = sim.bodies().vertices(Pose::to_poses(x, 3));

    std::vector<Pose> poses(sim.bodies().num_bodies());
    for (size_t i = 0; i < sim.bodies().num_bodies(); ++i) {
        poses[i] = Pose(x.segment<6>(6 * i));
    }

    sim.initialize_step();
    sim.update_collisions(x);
    REQUIRE(sim.normal_collisions().size() > 0);

    Eigen::VectorXd tmp = sim.barrier_potential().gradient(
        sim.normal_collisions(), sim.bodies(), V);

    Eigen::VectorXd fd_tmp;
    fd::finite_gradient(
        V.reshaped<Eigen::RowMajor>(),
        [&](const Eigen::VectorXd& x) -> double {
            return sim.barrier_potential()(
                sim.normal_collisions(), sim.bodies(),
                x.reshaped<Eigen::RowMajor>(V.rows(), V.cols()));
        },
        fd_tmp);

    CHECK(fd::compare_gradient(tmp, fd_tmp));
    if (!fd::compare_gradient(tmp, fd_tmp)) {
        std::cout << "analytic:\n" << tmp.transpose() << "\n\n";
        std::cout << "numerical:\n" << fd_tmp.transpose() << "\n\n";
    }

    // --- Check total gradient and hessian ------------------------------------

    Eigen::VectorXd g;
    sim.gradient(x, g);

    Eigen::VectorXd fd_g;
    fd::finite_gradient(
        x, [&](const Eigen::VectorXd& x) { return sim.value(x); }, fd_g);

    CHECK(fd::compare_gradient(g, fd_g));
    if (!fd::compare_gradient(g, fd_g)) {
        std::cout << "analytic:\n" << g << "\n\n";
        std::cout << "numerical:\n" << fd_g << "\n\n";
    }

    sim.set_project_to_psd(false); // FD comparison needs the exact Hessian
    Eigen::SparseMatrix<double> H_sparse;
    sim.hessian(x, H_sparse);
    const Eigen::MatrixXd H = H_sparse;

    Eigen::MatrixXd fd_H;
    fd::finite_jacobian(
        x,
        [&](const Eigen::VectorXd& x) {
            Eigen::VectorXd g_;
            sim.gradient(x, g_);
            return g_;
        },
        fd_H);

    CHECK(fd::compare_hessian(H, fd_H));
    if (!fd::compare_hessian(H, fd_H)) {
        std::cout << "analytic:\n" << H << "\n\n";
        std::cout << "numerical:\n" << fd_H << "\n\n";
    }

    sim.reset();
}
} // namespace ipc::tests

TEST_CASE("Rigid body simulator", "[.][rigid]")
{
    Eigen::MatrixXd V_bunny;
    Eigen::MatrixXi E_bunny, F_bunny;
    REQUIRE(tests::load_mesh("bunny (lowpoly).ply", V_bunny, E_bunny, F_bunny));

    Eigen::MatrixXd V_bowl;
    Eigen::MatrixXi E_bowl, F_bowl;
    REQUIRE(tests::load_mesh("bowl.ply", V_bowl, E_bowl, F_bowl));

    std::vector<Pose> initial_poses(3);
    initial_poses[0].position = Eigen::Vector3d(1.0, 1.5, 0);
    initial_poses[0].rotation = Eigen::Vector3d::Zero();
    initial_poses[1].position = Eigen::Vector3d(-1.0, 2.0, 0.0);
    initial_poses[1].rotation = Eigen::Vector3d(0.0, igl::PI / 4, 0.0);
    initial_poses[2].position = Eigen::Vector3d(0.0, 1.1, 0.0);
    initial_poses[2].rotation = Eigen::Vector3d(0.0, 0.0, 0.0);

    auto bodies = RigidBodies::build_from_meshes(
        std::vector<Eigen::MatrixXd> { V_bunny, V_bunny, V_bowl },
        std::vector<Eigen::MatrixXi> { E_bunny, E_bunny, E_bowl },
        std::vector<Eigen::MatrixXi> { F_bunny, F_bunny, F_bowl },
        /*densisties=*/ { { 1000.0, 1000.0, 1000.0 } }, initial_poses);
    bodies->planes.emplace_back(Eigen::Vector3d(0, 1, 0), 0);

    double dt = 0.1;
    double tend = 10.0;
    int n_steps = int(tend / dt);
    Simulator sim(bodies, initial_poses, dt);

    // --- Check simulator gradient and hessian --------------------------------
    SECTION("Simulator gradient and hessian")
    {
        ipc::tests::simulator_test(sim);
    }

    // --- Test stepping and running the simulation ----------------------------

    SECTION("Simulator stepping and running")
    {
        int n_calls = 0;
        auto callback = [&](bool success) {
            n_calls++;
            CHECK(success); // Ensure all steps succeed
        };

        REQUIRE(sim.step());
        CHECK(sim.t() == Catch::Approx(dt));
        REQUIRE(sim.run(tend, callback));
        CHECK(n_calls == n_steps - 1);
        CHECK(sim.t() == Catch::Approx(tend));
        CHECK(!sim.run(
            tend)); // Simulation already completed, should return false
        CHECK(n_calls == n_steps - 1); // Callback should not be called again

#ifdef IPC_TOOLKIT_WITH_GLTF
        write_gltf("simulator_test.glb", *bodies, sim.rigid_pose_history(), dt);
#endif

        sim.reset();
        CHECK(sim.t() == 0.0);
    }
}

TEST_CASE("Rigid simulator ballistic motion", "[rigid][simulator]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    tests::load_mesh("cube.ply", V, E, F);

    std::vector<Pose> initial_poses(1, Pose::Identity(3));
    initial_poses[0].position = Eigen::Vector3d(0, 100.0, 0); // free flight

    auto bodies = RigidBodies::build_from_meshes(
        { V }, { E }, { F }, { 1000.0 }, initial_poses);

    std::vector<Pose> initial_velocities(1, Pose::Identity(3));
    initial_velocities[0].position = Eigen::Vector3d(1.0, 2.0, -0.5); // v₀
    initial_velocities[0].rotation = Eigen::Vector3d(0, 0, 0.7);      // ω

    Simulator::Settings settings;
    settings.solver_params["grad_norm_tol"] = 1e-8;

    const double dt = 0.01;
    Simulator sim(bodies, initial_poses, initial_velocities, dt, settings);

    // NOTE: build_from_meshes folds R₀ into initial_poses, so capture the
    // world rotation after construction.
    const Eigen::Vector3d p0 = initial_poses[0].position;
    const Eigen::Matrix3d Q0 = initial_poses[0].rotation_matrix();
    const Eigen::Vector3d v0 = initial_velocities[0].position;
    const Eigen::Vector3d omega = initial_velocities[0].rotation;
    const Eigen::Vector3d g = settings.gravity;

    REQUIRE(sim.step());
    const AffinePose pose_1 = sim.pose_history().back()[0];

    // Implicit Euler translation: p₁ = p₀ + dt v₀ + dt² g
    const Eigen::Vector3d p1_expected = p0 + dt * v0 + dt * dt * g;
    CHECK(
        (pose_1.position - p1_expected).norm()
        == Catch::Approx(0).margin(1e-6));

    // Rotation: for a cube J ∝ I, so the minimizer of the inertial term is the
    // polar factor of Q̂ = (I + dt[ω]×) Q₀, i.e., a rotation about ω̂ by
    // atan(dt‖ω‖).
    const Eigen::Matrix3d Q1_expected =
        Eigen::AngleAxisd(std::atan(dt * omega.norm()), omega.normalized())
            .toRotationMatrix()
        * Q0;
    CHECK(
        (pose_1.rotation - Q1_expected).norm()
        == Catch::Approx(0).margin(1e-5));

    // rigid_pose_history() is the exact log map in rigid mode
    const Pose rigid_pose_1 = sim.rigid_pose_history().back()[0];
    CHECK(rigid_pose_1.position == pose_1.position);
    CHECK(
        (rigid_pose_1.rotation_matrix() - pose_1.rotation).norm()
        == Catch::Approx(0).margin(1e-12));

    // Multi-step translation closed form:
    // pₙ = p₀ + n dt v₀ + dt² g n(n+1)/2
    const int n = 10;
    for (int i = 1; i < n; i++) {
        REQUIRE(sim.step());
    }
    const Eigen::Vector3d pn_expected =
        p0 + n * dt * v0 + dt * dt * g * (n * (n + 1) / 2.0);
    CHECK(
        (sim.pose_history().back()[0].position - pn_expected).norm()
        == Catch::Approx(0).margin(1e-5));

    // reset() must restore the initial velocities (and rebuild the body
    // forces with the new time integrator).
    sim.reset();
    CHECK(sim.t() == 0.0);
    REQUIRE(sim.step());
    const AffinePose pose_1_again = sim.pose_history().back()[0];
    CHECK(
        (pose_1_again.position - pose_1.position).norm()
        == Catch::Approx(0).margin(1e-10));
    CHECK(
        (pose_1_again.rotation - pose_1.rotation).norm()
        == Catch::Approx(0).margin(1e-10));
}

TEST_CASE("Rigid simulator settings", "[rigid][simulator]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    tests::load_mesh("cube.ply", V, E, F);

    std::vector<Pose> initial_poses(1, Pose::Identity(3));
    auto bodies = RigidBodies::build_from_meshes(
        { V }, { E }, { F }, { 1000.0 }, initial_poses);

    Simulator::Settings settings;
    settings.dhat = 0.05;
    settings.gravity = Eigen::Vector3d(0, -1.62, 0); // moon gravity
    settings.solver_params["max_iterations"] = 42;

    Simulator sim(bodies, initial_poses, /*dt=*/0.01, settings);

    CHECK(sim.barrier_potential().dhat() == 0.05);
    CHECK(sim.gravity() == Eigen::Vector3d(0, -1.62, 0));
    CHECK(sim.settings().solver_params["max_iterations"] == 42);

    sim.set_gravity(Eigen::Vector3d(0, -9.81, 0));
    CHECK(sim.gravity() == Eigen::Vector3d(0, -9.81, 0));
}

// Regression test: a bunny dropped off-center onto the bowl rim used to tunnel
// through the thin shell because the rigid continuous broad phase missed the
// bunny-bowl candidate pair. The narrow phase (CCD) is only as good as the
// candidates it is given, so the simulation must stay intersection-free.
TEST_CASE("Rigid simulator no tunneling into bowl", "[rigid][simulator]")
{
    Eigen::MatrixXd V_bunny, V_bowl;
    Eigen::MatrixXi E_bunny, F_bunny, E_bowl, F_bowl;
    REQUIRE(tests::load_mesh("bunny (lowpoly).ply", V_bunny, E_bunny, F_bunny));
    REQUIRE(tests::load_mesh("bowl.ply", V_bowl, E_bowl, F_bowl));

    // Bunny off to the side so it falls onto the sloped rim/wall of the bowl.
    std::vector<Pose> initial_poses(2, Pose::Identity(3));
    initial_poses[0].position = Eigen::Vector3d(1.0, 1.5, 0.0); // bunny
    initial_poses[1].position = Eigen::Vector3d(0.0, 1.1, 0.0); // bowl

    auto bodies = RigidBodies::build_from_meshes(
        { V_bunny, V_bowl }, { E_bunny, E_bowl }, { F_bunny, F_bowl },
        { 1000.0, 1000.0 }, initial_poses);
    bodies->planes.emplace_back(Eigen::Vector3d(0, 1, 0), 0); // ground y = 0

    Simulator sim(bodies, initial_poses, /*dt=*/1.0 / 60.0);

    // First contact with the bowl occurs around t ≈ 0.5; run past it.
    for (int i = 0; i < 60; ++i) {
        REQUIRE(sim.step());
        const Eigen::MatrixXd V =
            bodies->vertices(sim.rigid_pose_history().back());
        CHECK(!has_intersections(*bodies, V));
    }
}

TEST_CASE("Rigid simulator BDF order", "[rigid][simulator]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    tests::load_mesh("cube.ply", V, E, F);

    std::vector<Pose> initial_poses(1, Pose::Identity(3));
    initial_poses[0].position = Eigen::Vector3d(0, 100.0, 0);

    auto bodies = RigidBodies::build_from_meshes(
        { V }, { E }, { F }, { 1000.0 }, initial_poses);

    std::vector<Pose> initial_velocities(1, Pose::Identity(3));
    initial_velocities[0].position = Eigen::Vector3d(1.0, -2.0, 0.5); // v₀

    // Constant-velocity free flight (no gravity, no torque) is a linear
    // trajectory, which any BDF order integrates exactly.
    const int order = GENERATE(1, 2, 3, 4);
    Simulator::Settings settings;
    settings.bdf_order = order;
    settings.gravity.setZero();
    settings.solver_params["grad_norm_tol"] = 1e-12;

    const double dt = 0.01;
    const Eigen::Vector3d p0 = initial_poses[0].position;
    const Eigen::Vector3d v0 = initial_velocities[0].position;

    Simulator sim(bodies, initial_poses, initial_velocities, dt, settings);

    const int n = 10;
    for (int i = 0; i < n; ++i) {
        REQUIRE(sim.step());
    }

    // pₙ = p₀ + n·dt·v₀ exactly, for any BDF order.
    const Eigen::Vector3d pn = sim.pose_history().back()[0].position;
    CHECK((pn - (p0 + n * dt * v0)).norm() == Catch::Approx(0).margin(1e-9));
}
