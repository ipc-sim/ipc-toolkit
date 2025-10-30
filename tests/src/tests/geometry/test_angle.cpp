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

    double expected_angle;
    SECTION("Various angles")
    {
        const double rot = deg2rad(
            GENERATE(1, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 179));
        CAPTURE(rot);

        Eigen::AngleAxisd rotation(
            (igl::PI - rot) / 2.0, Eigen::Vector3d::UnitY());
        x2 = rotation * x2;
        x3 = rotation.inverse() * x3;

        expected_angle = rot - igl::PI;
    }

    SECTION("Case 1")
    {
        x0 << -0.015247385606936873, 1.1187166216183693, -0.09508569727171673;
        x1 << -0.017971426013627917, 1.1229485012592226, -0.0934495693604115;
        x2 << -0.021023579437439658, 1.1190719774332136, -0.09490934871219975;
        x3 << -0.014473605843359506, 1.1216871870663812, -0.09395608203835042;
        expected_angle = 0;
    }

    double angle = dihedral_angle(x0, x1, x2, x3);
    CHECK(angle == Catch::Approx(expected_angle).margin(1e-8));

    // --- Check the gradient against finite differences ---

    Eigen::VectorXd x(12);
    x << x0, x1, x2, x3;

    Eigen::VectorXd grad = dihedral_angle_gradient(x0, x1, x2, x3);
    Eigen::VectorXd fd_grad(12);
    fd::finite_gradient(
        x,
        [](const Eigen::VectorXd& x_fd) -> double {
            double theta = dihedral_angle(
                x_fd.segment<3>(0), x_fd.segment<3>(3), x_fd.segment<3>(6),
                x_fd.segment<3>(9));
            return theta;
        },
        fd_grad);

    REQUIRE(grad.array().isFinite().all());
    CHECK(fd::compare_gradient(grad, fd_grad));
    if (!fd::compare_gradient(grad, fd_grad)) {
        std::cout << "   Gradient:\n" << grad.transpose() << std::endl;
        std::cout << "FD Gradient:\n" << fd_grad.transpose() << std::endl;
    }
}