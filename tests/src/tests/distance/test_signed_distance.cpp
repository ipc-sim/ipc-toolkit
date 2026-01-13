#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/catch_approx.hpp>

#include <ipc/distance/signed/point_line.hpp>
#include <ipc/distance/signed/point_plane.hpp>
#include <ipc/distance/signed/line_line.hpp>

#include <iostream>
#include <finitediff.hpp>

using namespace ipc;

TEST_CASE(
    "Signed distance point-plane derivatives", "[normal][gradient][hessian]")
{
    Eigen::Vector3d p(0.25, 0.25, 1.0);
    Eigen::Vector3d a(0, 0, 0);
    Eigen::Vector3d b(1, 0, 0);
    Eigen::Vector3d c(0, 1, 0);
    Eigen::VectorXd x(12);

    const int case_i = GENERATE(range(0, 10));
    if (case_i > 0) {
        p.setRandom();
        a.setRandom();
        b.setRandom();
        c.setRandom();
    }

    x << p, a, b, c;

    // Check gradient using finite differences
    Vector12d gradient = point_plane_signed_distance_gradient(p, a, b, c);
    Eigen::VectorXd fd_gradient;
    fd::finite_gradient(
        x,
        [](const Eigen::VectorXd& x_fd) -> double {
            return point_plane_signed_distance(
                x_fd.segment<3>(0), x_fd.segment<3>(3), x_fd.segment<3>(6),
                x_fd.segment<3>(9));
        },
        fd_gradient);
    CHECK(fd::compare_gradient(gradient, fd_gradient, 1e-6));
    if (!fd::compare_gradient(gradient, fd_gradient, 1e-6)) {
        std::cout << "Gradient:\n" << gradient.transpose() << std::endl;
        std::cout << "FD Gradient:\n" << fd_gradient.transpose() << std::endl;
    }

    // Check hessian using finite differences
    Matrix12d hessian = point_plane_signed_distance_hessian(p, a, b, c);
    Eigen::MatrixXd fd_hessian;
    fd::finite_jacobian(
        x,
        [](const Eigen::VectorXd& x_fd) -> Eigen::VectorXd {
            return point_plane_signed_distance_gradient(
                x_fd.segment<3>(0), x_fd.segment<3>(3), x_fd.segment<3>(6),
                x_fd.segment<3>(9));
        },
        fd_hessian);
    CHECK(fd::compare_jacobian(hessian, fd_hessian, 1e-6));
    if (!fd::compare_jacobian(hessian, fd_hessian, 1e-6)) {
        std::cout << "Hessian:\n" << hessian << std::endl;
        std::cout << "FD Hessian:\n" << fd_hessian << std::endl;
    }
}

TEST_CASE(
    "Signed distance line-line derivatives", "[normal][gradient][hessian]")
{
    Eigen::Vector3d ea0(0.25, 0.25, 1.0);
    Eigen::Vector3d ea1(0, 0, 0);
    Eigen::Vector3d eb0(1, 0, 0);
    Eigen::Vector3d eb1(0, 1, 0);
    Eigen::VectorXd x(12);

    const int case_i = GENERATE(range(0, 10));
    while (case_i > 0 && line_line_signed_distance(ea0, ea1, eb0, eb1) == 0.0) {
        ea0.setRandom();
        ea1.setRandom();
        eb0.setRandom();
        eb1.setRandom();
    }

    x << ea0, ea1, eb0, eb1;

    // Check gradient using finite differences
    Vector12d gradient = line_line_signed_distance_gradient(ea0, ea1, eb0, eb1);
    Eigen::VectorXd fd_gradient;
    fd::finite_gradient(
        x,
        [](const Eigen::VectorXd& x_fd) -> double {
            return line_line_signed_distance(
                x_fd.segment<3>(0), x_fd.segment<3>(3), x_fd.segment<3>(6),
                x_fd.segment<3>(9));
        },
        fd_gradient);
    CHECK(fd::compare_gradient(gradient, fd_gradient, 1e-6));
    if (!fd::compare_gradient(gradient, fd_gradient, 1e-6)) {
        std::cout << "Gradient:\n" << gradient.transpose() << std::endl;
        std::cout << "FD Gradient:\n" << fd_gradient.transpose() << std::endl;
    }

    // Check hessian using finite differences
    Matrix12d hessian = line_line_signed_distance_hessian(ea0, ea1, eb0, eb1);
    Eigen::MatrixXd fd_hessian;
    fd::finite_jacobian(
        x,
        [](const Eigen::VectorXd& x_fd) -> Eigen::VectorXd {
            return line_line_signed_distance_gradient(
                x_fd.segment<3>(0), x_fd.segment<3>(3), x_fd.segment<3>(6),
                x_fd.segment<3>(9));
        },
        fd_hessian);
    CHECK(fd::compare_jacobian(hessian, fd_hessian, 1e-6));
    if (!fd::compare_jacobian(hessian, fd_hessian, 1e-6)) {
        std::cout << "Hessian:\n" << hessian << std::endl;
        std::cout << "FD Hessian:\n" << fd_hessian << std::endl;
    }
}

TEST_CASE(
    "Signed distance point-line derivatives", "[normal][gradient][hessian]")
{
    Eigen::Vector2d p(0.25, 0.25);
    Eigen::Vector2d a(0, 0);
    Eigen::Vector2d b(1, 0);
    Eigen::VectorXd x(6);

    const int case_i = GENERATE(range(0, 10));
    if (case_i > 0) {
        p.setRandom();
        a.setRandom();
        b.setRandom();
    }

    x << p, a, b;

    // Check gradient using finite differences
    Vector6d gradient = point_line_signed_distance_gradient(p, a, b);
    Eigen::VectorXd fd_gradient;
    fd::finite_gradient(
        x,
        [](const Eigen::VectorXd& x_fd) -> double {
            return point_line_signed_distance(
                x_fd.segment<2>(0), x_fd.segment<2>(2), x_fd.segment<2>(4));
        },
        fd_gradient);
    CHECK(fd::compare_gradient(gradient, fd_gradient, 1e-6));
    if (!fd::compare_gradient(gradient, fd_gradient, 1e-6)) {
        std::cout << "Gradient:\n" << gradient.transpose() << std::endl;
        std::cout << "FD Gradient:\n" << fd_gradient.transpose() << std::endl;
    }

    // Check hessian using finite differences
    Matrix6d hessian = point_line_signed_distance_hessian(p, a, b);
    Eigen::MatrixXd fd_hessian;
    fd::finite_jacobian(
        x,
        [](const Eigen::VectorXd& x_fd) -> Eigen::VectorXd {
            return point_line_signed_distance_gradient(
                x_fd.segment<2>(0), x_fd.segment<2>(2), x_fd.segment<2>(4));
        },
        fd_hessian);
    CHECK(fd::compare_jacobian(hessian, fd_hessian, 1e-6));
    if (!fd::compare_jacobian(hessian, fd_hessian, 1e-6)) {
        std::cout << "Hessian:\n" << hessian << std::endl;
        std::cout << "FD Hessian:\n" << fd_hessian << std::endl;
    }
}