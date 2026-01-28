// Test the mass utilities.

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <tests/utils.hpp>

#include <ipc/dynamics/rigid/rigid_body.hpp>

using namespace ipc;
using namespace ipc::rigid;

TEST_CASE("Rigid body construction", "[rigid]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    REQUIRE(tests::load_mesh("cube.ply", V, E, F));

    const double L = 0.5, W = 1.0, H = 2.0;
    V.col(0).array() *= L;
    V.col(1).array() *= W;
    V.col(2).array() *= H;

    const double density = GENERATE(1.0, 2.0, 3.0);

    Pose input_pose = Pose::Zero(3);
    Pose initial_pose = Pose::Zero(3);

    SECTION("No modification") { }
    SECTION("Input pose")
    {
        input_pose.position = Eigen::Vector3d(1.0, 2.0, 3.0);
        input_pose.rotation = Eigen::Vector3d(0.1, 0.2, 0.3);
    }
    SECTION("Initial pose")
    {
        initial_pose.position = Eigen::Vector3d(1.0, 2.0, 3.0);
        initial_pose.rotation = Eigen::Vector3d(0.1, 0.2, 0.3);
    }
    SECTION("Input and initial pose")
    {
        input_pose.position = Eigen::Vector3d(1.0, 2.0, 3.0);
        input_pose.rotation = Eigen::Vector3d(0.1, 0.2, 0.3);
        initial_pose.position = Eigen::Vector3d(4.0, 5.0, 6.0);
        initial_pose.rotation = Eigen::Vector3d(0.4, 0.5, 0.6);
    }

    Eigen::MatrixXd modified_V = input_pose.transform_vertices(V);

    const double m = density * L * W * H; // unit mass per voxel
    const double Ixx = m * (W * W + H * H) / 12.0;
    const double Iyy = m * (L * L + H * H) / 12.0;
    const double Izz = m * (L * L + W * W) / 12.0;
    Eigen::Vector3d I(Ixx, Iyy, Izz);
    // Sort:
    if (I(0) > I(1)) {
        std::swap(I(0), I(1));
    }
    if (I(1) > I(2)) {
        std::swap(I(1), I(2));
    }
    if (I(0) > I(1)) {
        std::swap(I(0), I(1));
    }

    Pose modified_pose = initial_pose;
    const RigidBody body(modified_V, E, F, density, modified_pose);
    REQUIRE(modified_V.rows() == V.rows());
    REQUIRE(modified_V.cols() == V.cols());

    {
        CAPTURE(body.moment_of_inertia().transpose(), I.transpose());
        REQUIRE(body.mass() == Catch::Approx(m).margin(1e-8));
        REQUIRE(body.moment_of_inertia().isApprox(I));
    }

    // The modified position should be the initial position plus the input
    // position
    {
        CAPTURE(
            input_pose.position.transpose(), initial_pose.position.transpose(),
            modified_pose.position.transpose());
        REQUIRE(modified_pose.position.isApprox(
            input_pose.position + initial_pose.position));
    }

    // The modified rotation should be the initial rotation times the input
    // rotation
    {
        const Eigen::MatrixXd _V = modified_pose.transform_vertices(modified_V);

        const Eigen::MatrixXd expected_V =
            (initial_pose * input_pose).transform_vertices(V);

        CHECK(_V.isApprox(expected_V, 1e-8));
    }
}