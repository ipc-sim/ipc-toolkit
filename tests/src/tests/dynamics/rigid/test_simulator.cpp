// Test the mass utilities.

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <tests/utils.hpp>

#include <ipc/dynamics/rigid/simulator.hpp>
#include <ipc/dynamics/rigid/rigid_bodies.hpp>
#include <ipc/io/write_gltf.hpp>
#include <ipc/potentials/barrier_potential.hpp>

#include <finitediff.hpp>
#include <igl/PI.h>

#include <iostream>

using namespace ipc;
using namespace ipc::rigid;

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

    const Eigen::VectorXd g = sim.gradient(x);

    Eigen::VectorXd fd_g;
    fd::finite_gradient(
        x, [&](const Eigen::VectorXd& x) { return sim.energy(x); }, fd_g);

    CHECK(fd::compare_gradient(g, fd_g));
    if (!fd::compare_gradient(g, fd_g)) {
        std::cout << "analytic:\n" << g << "\n\n";
        std::cout << "numerical:\n" << fd_g << "\n\n";
    }

    const Eigen::MatrixXd H = sim.hessian(x, false);

    Eigen::MatrixXd fd_H;
    fd::finite_jacobian(
        x, [&](const Eigen::VectorXd& x) { return sim.gradient(x); }, fd_H);

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

        write_gltf("simulator_test.glb", *bodies, sim.pose_history(), dt);

        sim.reset();
        CHECK(sim.t() == 0.0);
    }
}