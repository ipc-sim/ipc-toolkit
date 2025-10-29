#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/catch_approx.hpp>

#include <ipc/geometry/angle.hpp>

#include <igl/PI.h>
#include <finitediff.hpp>

#include <iostream>

using namespace ipc;

constexpr double deg2rad(double degrees) { return degrees * igl::PI / 180.0; }

TEST_CASE("Dihedral angle and gradient", "[angle][dihedral]")
{
    Eigen::Vector3d x0(0, -1, 0);
    Eigen::Vector3d x1(0, 1, 0);
    Eigen::Vector3d x2(0.5, 0, 0);
    Eigen::Vector3d x3(-0.5, 0, 0);

    double expected_angle = deg2rad(GENERATE(0, 15, 30, 45, 60, 75, 90));

    Eigen::AngleAxisd rotation(
        (igl::PI - expected_angle) / 2.0, Eigen::Vector3d::UnitY());
    x2 = rotation * x2;
    x3 = rotation.inverse() * x3;

    double angle = dihedral_angle(x0, x1, x2, x3);

    if (expected_angle < igl::PI / 2) {
        expected_angle = igl::PI - expected_angle;
    }

    CHECK(angle == Catch::Approx(expected_angle));

    // --- Check the gradient against finite differences ---

    Eigen::VectorXd x(12);
    x << x0, x1, x2, x3;

    Eigen::VectorXd grad = dihedral_angle_gradient(x0, x1, x2, x3);
    Eigen::VectorXd fd_grad(12);
    fd::finite_gradient(
        x,
        [](const Eigen::VectorXd& x_fd) -> double {
            return dihedral_angle(
                x_fd.segment<3>(0), x_fd.segment<3>(3), x_fd.segment<3>(6),
                x_fd.segment<3>(9));
        },
        fd_grad);

    CHECK(fd::compare_gradient(grad, fd_grad));
    if (!fd::compare_gradient(grad, fd_grad)) {
        std::cout << "   Gradient:\n" << grad.transpose() << std::endl;
        std::cout << "FD Gradient:\n" << fd_grad.transpose() << std::endl;
    }
}