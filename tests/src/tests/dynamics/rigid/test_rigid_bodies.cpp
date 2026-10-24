// Test the mass utilities.

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <finitediff.hpp>
#include <tests/utils.hpp>
#include <iostream>

#include <ipc/dynamics/rigid/rigid_bodies.hpp>

using namespace ipc;
using namespace ipc::rigid;

TEST_CASE("RigidBodies::to_rigid_dof", "[rigid]")
{
    const int num_vertices_per_body = 10;

    std::vector<Pose> initial_poses(2, Pose::Identity(3));
    auto bodies = RigidBodies::build_from_meshes(
        std::vector<Eigen::MatrixXd> {
            Eigen::MatrixXd::Random(num_vertices_per_body, 3),
            Eigen::MatrixXd::Random(num_vertices_per_body, 3),
        },
        std::vector<Eigen::MatrixXi>(2), std::vector<Eigen::MatrixXi>(2),
        std::vector<double>(2, 1.0), initial_poses);

    auto f = [&](const Eigen::VectorXd& x) -> double {
        assert(x.size() == 12);
        auto Va = bodies->body_vertices(0, Pose(x.head<6>()));
        auto Vb = bodies->body_vertices(1, Pose(x.tail<6>()));
        return (Va.transpose() * Vb).trace();
    };

    Eigen::VectorXd x = Eigen::VectorXd::Random(12);

    // --- Gradient ------------------------------------------------------------

    Eigen::VectorXd df_dV(bodies->ndof());
    df_dV.head(3 * num_vertices_per_body) =
        bodies->body_vertices(1, Pose(x.tail<6>())).reshaped<Eigen::RowMajor>();
    df_dV.tail(3 * num_vertices_per_body) =
        bodies->body_vertices(0, Pose(x.head<6>())).reshaped<Eigen::RowMajor>();

    Eigen::VectorXd df_dx =
        bodies->to_rigid_dof({ Pose(x.head<6>()), Pose(x.tail<6>()) }, df_dV);

    Eigen::VectorXd df_dx_fd;
    fd::finite_gradient(x, f, df_dx_fd);

    CHECK(fd::compare_gradient(df_dx, df_dx_fd));
    if (!fd::compare_gradient(df_dx, df_dx_fd)) {
        std::cout << "analytic:\n" << df_dx << "\n\n";
        std::cout << "numerical:\n" << df_dx_fd << "\n\n";
    }

    // --- Hessian -------------------------------------------------------------

    Eigen::MatrixXd d2f_dV2 =
        Eigen::MatrixXd::Zero(bodies->ndof(), bodies->ndof());
    d2f_dV2.topRightCorner(3 * num_vertices_per_body, 3 * num_vertices_per_body)
        .setIdentity();
    d2f_dV2
        .bottomLeftCorner(3 * num_vertices_per_body, 3 * num_vertices_per_body)
        .setIdentity();

    Eigen::MatrixXd d2f_dx2 = bodies->to_rigid_dof(
        { Pose(x.head<6>()), Pose(x.tail<6>()) }, df_dV, d2f_dV2.sparseView());

    Eigen::MatrixXd d2f_dx2_fd;
    fd::finite_hessian(x, f, d2f_dx2_fd);

    CHECK(fd::compare_hessian(d2f_dx2, d2f_dx2_fd));
    if (!fd::compare_hessian(d2f_dx2, d2f_dx2_fd)) {
        std::cout << "analytic:\n" << d2f_dx2 << "\n\n";
        std::cout << "numerical:\n" << d2f_dx2_fd << "\n\n";
    }
}
TEST_CASE(
    "Converting a plane keeps the other bodies' properties aligned", "[rigid]")
{
    Eigen::MatrixXd V_cube;
    Eigen::MatrixXi E_cube, F_cube;
    REQUIRE(tests::load_mesh("cube.ply", V_cube, E_cube, F_cube));

    // A 2×2 square in the y = 0 plane, listed before the cubes.
    Eigen::MatrixXd V_plane(4, 3);
    V_plane << -1, 0, -1, 1, 0, -1, 1, 0, 1, -1, 0, 1;
    Eigen::MatrixXi F_plane(2, 3);
    F_plane << 0, 2, 1, 0, 3, 2;
    Eigen::MatrixXi E_plane(5, 2);
    E_plane << 0, 1, 1, 2, 2, 0, 2, 3, 3, 0;

    VectorMax6b slider(6);
    slider << false, true, true, true, true, true;

    std::vector<Pose> initial_poses(3, Pose::Identity(3));
    initial_poses[1].position = Eigen::Vector3d(0, 1, 0);
    initial_poses[2].position = Eigen::Vector3d(3, 1, 0);

    const auto bodies = RigidBodies::build_from_meshes(
        { V_plane, V_cube, V_cube }, { E_plane, E_cube, E_cube },
        { F_plane, F_cube, F_cube }, { 1.0, 500.0, 2000.0 }, initial_poses,
        /*convert_planes=*/true, { VectorMax6b(), VectorMax6b(), slider },
        { "ground", "light", "heavy" });

    REQUIRE(bodies->planes.size() == 1);
    REQUIRE(bodies->num_bodies() == 2);
    REQUIRE(initial_poses.size() == 2);

    // Each cube keeps the density, name, DOF mask, and pose it was given.
    CHECK((*bodies)[0].density() == Catch::Approx(500.0));
    CHECK((*bodies)[1].density() == Catch::Approx(2000.0));
    CHECK(bodies->body_name(0) == "light");
    CHECK(bodies->body_name(1) == "heavy");
    CHECK(!(*bodies)[0].is_dof_fixed().any());
    CHECK((*bodies)[1].is_dof_fixed() == slider);
    CHECK(initial_poses[0].position.isApprox(Eigen::Vector3d(0, 1, 0)));
    CHECK(initial_poses[1].position.isApprox(Eigen::Vector3d(3, 1, 0)));
}
