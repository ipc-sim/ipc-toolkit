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