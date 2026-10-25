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

    Pose input_pose = Pose::identity(3);
    Pose initial_pose = Pose::identity(3);

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
            (initial_pose * input_pose).position));
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
TEST_CASE(
    "Changing a body's density matches building it at that density", "[rigid]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    REQUIRE(tests::load_mesh("cube.ply", V, E, F));
    V.col(0).array() *= 2.0; // anisotropic so the inertia is nontrivial

    const double density = GENERATE(0.5, 7850.0);

    Eigen::MatrixXd V_default = V, V_expected = V;
    Pose pose_default = Pose::identity(3), pose_expected = Pose::identity(3);
    RigidBody body(V_default, E, F, /*density=*/1000.0, pose_default);
    const RigidBody expected(V_expected, E, F, density, pose_expected);

    body.set_density(density);

    CHECK(body.density() == Catch::Approx(density));
    CHECK(body.volume() == Catch::Approx(expected.volume()));
    CHECK(body.mass() == Catch::Approx(expected.mass()));
    CHECK(body.moment_of_inertia().isApprox(expected.moment_of_inertia()));
    CHECK(body.second_moment().diagonal().isApprox(
        expected.second_moment().diagonal()));
    // The geometry and principal frame do not depend on density.
    CHECK(body.principal_rotation().isApprox(expected.principal_rotation()));
    CHECK(pose_default.rotation.isApprox(pose_expected.rotation));
}

TEST_CASE(
    "A body spinning about one world axis has that axis's inertia", "[rigid]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi E, F;
    REQUIRE(tests::load_mesh("cube.ply", V, E, F));
    const double L = 0.5, W = 1.0, H = 2.0;
    V.col(0).array() *= L;
    V.col(1).array() *= W;
    V.col(2).array() *= H;

    const int free_axis = GENERATE(0, 1, 2);
    CAPTURE(free_axis);

    VectorMax6b is_dof_fixed = VectorMax6b::Zero(6);
    is_dof_fixed.tail<3>().setOnes();
    is_dof_fixed(3 + free_axis) = false;

    const double density = 1000.0;
    Pose pose = Pose::identity(3);
    const RigidBody body(V, E, F, density, pose, is_dof_fixed);

    // The second moment ∫ρ x̄x̄ᵀ of the box, whose two entries orthogonal to
    // the free axis sum to the moment of inertia about that axis.
    const double m = density * L * W * H;
    const Eigen::Vector3d second_moment =
        m / 12.0 * Eigen::Vector3d(L * L, W * W, H * H);

    REQUIRE(body.second_moment().diagonal().size() == 3);
    CHECK(body.second_moment().diagonal().isApprox(second_moment));
    const double I_free = second_moment.sum() - second_moment(free_axis);
    CHECK(
        body.moment_of_inertia()(free_axis)
        == Catch::Approx(I_free).epsilon(1e-8));
}
