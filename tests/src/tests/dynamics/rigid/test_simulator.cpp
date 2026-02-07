// Test the mass utilities.

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <tests/utils.hpp>

#include <ipc/dynamics/rigid/simulator.hpp>
#include <ipc/dynamics/rigid/rigid_bodies.hpp>

#include <igl/PI.h>
#include <memory>

using namespace ipc;
using namespace ipc::rigid;

TEST_CASE("Rigid body simulator", "[rigid]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    REQUIRE(tests::load_mesh("bunny (lowpoly).ply", V, E, F));

    std::vector<Pose> initial_poses(2);
    initial_poses[0].position = Eigen::Vector3d(0, 0.5, 0);
    initial_poses[0].rotation = Eigen::Vector3d::Random();
    initial_poses[1].position = Eigen::Vector3d(2.0, 0.5, 0.0);
    initial_poses[1].rotation = Eigen::Vector3d(0.0, igl::PI / 4, 0.0);

    Eigen::MatrixXd rest_positions(2 * V.rows(), 3);
    rest_positions.topRows(V.rows()) = V;
    rest_positions.bottomRows(V.rows()) = V;

    Eigen::MatrixXi edges(2 * E.rows(), 2);
    edges.topRows(E.rows()) = E;
    edges.bottomRows(E.rows()) = E.array() + V.rows();

    Eigen::MatrixXi faces(2 * F.rows(), 3);
    faces.topRows(F.rows()) = F;
    faces.bottomRows(F.rows()) = F.array() + V.rows();

    auto bodies = std::make_shared<RigidBodies>(RigidBodies::build_from_meshes(
        std::vector<Eigen::MatrixXd> { V, V },
        std::vector<Eigen::MatrixXi> { E, E },
        std::vector<Eigen::MatrixXi> { F, F },
        /*densisties=*/ { { 1.0, 1.0 } }, initial_poses));

    Simulator sim(bodies, initial_poses, 0.1);

    int n_calls = 0;
    auto callback = [&n_calls](bool success) {
        n_calls++;
        CHECK(success); // Ensure all steps succeed
    };

    CHECK(sim.step());
    CHECK(sim.t() == Catch::Approx(0.1));
    CHECK(sim.run(1.0, callback)); // t = 1.0
    CHECK(n_calls == 9);           // 9 steps from t=0.1 to t=1.0
    CHECK(sim.t() == Catch::Approx(1.0));
    CHECK(!sim.run(1.0)); // Simulation already completed, should return false
    CHECK(n_calls == 9);  // Callback should not be called again
    sim.reset();
    CHECK(sim.t() == 0.0);
}