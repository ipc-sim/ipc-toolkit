// The code that the "Rigid and Affine Body Dynamics" tutorial shows, so that we
// know it keeps working.

#include <Eigen/Core>
#include <catch2/catch_all.hpp>

#include <tests/utils.hpp>

#include <ipc/dynamics/rigid/rigid_bodies.hpp>
#include <simulator.hpp>

using namespace ipc;
using ipc::demo::Simulator;

TEST_CASE(
    "A dropped cube comes to rest on the ground",
    "[tutorial][simulator][rigid][affine]")
{
    const auto body_dynamics = GENERATE(
        Simulator::BodyDynamics::RIGID, Simulator::BodyDynamics::AFFINE);

    // A unit cube centered at the origin with vertices V, edges E, and faces F.
    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    REQUIRE(tests::load_mesh("cube.ply", V, E, F));

    // BEGIN_DROP_A_CUBE
    // Hold the cube 0.5 m above the ground (its center is at y = 1).
    std::vector<rigid::Pose> poses = { rigid::Pose(
        Eigen::Vector3d(0, 1, 0), Eigen::Vector3d::Zero()) };

    auto bodies = rigid::RigidBodies::build_from_meshes(
        { V }, { E }, { F }, /*densities=*/ { 1000.0 }, poses);

    // The ground is the plane y = 0.
    bodies->planes.emplace_back(
        Eigen::Vector3d::UnitY(), Eigen::Vector3d::Zero());

    Simulator::Settings settings;
    settings.body_dynamics = body_dynamics; // RIGID or AFFINE

    Simulator sim(bodies, poses, /*dt=*/0.01, settings);
    const bool success = sim.run(/*t_end=*/0.6);

    const std::vector<rigid::Pose> final_poses = sim.rigid_poses();
    // END_DROP_A_CUBE

    CHECK(success);

    // The cube has landed: it is on the ground (not through it), and it has
    // fallen at least most of the way down.
    const double height = final_poses[0].position.y();
    CHECK(height >= 0.5 - 1e-8);
    CHECK(height < 0.6);
}
