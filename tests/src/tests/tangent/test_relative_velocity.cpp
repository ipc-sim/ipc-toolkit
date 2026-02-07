#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <finitediff.hpp>
#include <iostream>

#include <ipc/tangent/closest_point.hpp>
#include <ipc/tangent/relative_velocity.hpp>

using namespace ipc;

TEST_CASE(
    "Point-triangle relative velocity",
    "[friction][point-triangle][relative_velocity]")
{
    Eigen::Vector3d p(0, 1, 0), t0(-1, 0, 1), t1(1, 0, 1), t2(0, 0, -1);
    Eigen::Vector3d dp = Eigen::Vector3d::Random(), dt0, dt1, dt2;
    dt0 = dt1 = dt2 = Eigen::Vector3d::Zero();

    Eigen::Vector2d barycentric_coords =
        point_triangle_closest_point(p, t0, t1, t2);

    Eigen::Vector3d relative_velocity =
        point_triangle_relative_velocity(dp, dt0, dt1, dt2, barycentric_coords);

    CHECK((dp - relative_velocity).norm() == Catch::Approx(0).margin(1e-12));
}

TEST_CASE(
    "Edge-edge relative velocity", "[friction][edge-edge][relative_velocity]")
{
    Eigen::Vector3d ea0(-1, 0, 0), ea1(1, 0, 0), eb0(0, 0, -1), eb1(0, 0, 1);
    Eigen::Vector3d dea0, dea1, deb0(0, 0, 0), deb1(0, 0, 0);
    dea0 = dea1 = Eigen::Vector3d::Random();

    Eigen::Vector2d barycentric_coords =
        edge_edge_closest_point(ea0, ea1, eb0, eb1);

    Eigen::Vector3d relative_velocity =
        edge_edge_relative_velocity(dea0, dea1, deb0, deb1, barycentric_coords);

    CHECK((dea0 - relative_velocity).norm() == Catch::Approx(0).margin(1e-12));
}

TEST_CASE(
    "Point-edge relative velocity", "[friction][point-edge][relative_velocity]")
{
    Eigen::Vector3d p(0, 1, 0), e0(-1, 0, 0), e1(1, 0, 0);
    Eigen::Vector3d dp = Eigen::Vector3d::Random(), de0(0, 0, 0), de1(0, 0, 0);

    double alpha = point_edge_closest_point(p, e0, e1);

    Eigen::Vector3d relative_velocity =
        point_edge_relative_velocity(dp, de0, de1, alpha);
    CHECK((dp - relative_velocity).norm() == Catch::Approx(0).margin(1e-12));
}

TEST_CASE(
    "Point-point relative velocity",
    "[friction][point-point][relative_velocity]")
{
    Eigen::Vector3d dp0 = Eigen::Vector3d::Random(), dp1(0, 0, 0);

    Eigen::Vector3d relative_velocity = point_point_relative_velocity(dp0, dp1);
    CHECK((dp0 - relative_velocity).norm() == Catch::Approx(0).margin(1e-12));
}

TEST_CASE(
    "Point-edge relative velocity in 2D",
    "[friction][point-edge][relative_velocity][2D]")
{
    Eigen::Vector2d p(0, 1), e0(-1, 0), e1(1, 0);
    Eigen::Vector2d dp = Eigen::Vector2d::Random(), de0(0, 0), de1(0, 0);

    double alpha = point_edge_closest_point(p, e0, e1);

    Eigen::Vector2d relative_velocity =
        point_edge_relative_velocity(dp, de0, de1, alpha);
    CHECK((dp - relative_velocity).norm() == Catch::Approx(0).margin(1e-12));
}

TEST_CASE(
    "Point-point relative velocity in 2D",
    "[friction][point-point][relative_velocity][2D]")
{
    Eigen::Vector2d dp0 = Eigen::Vector2d::Random(), dp1(0, 0);

    Eigen::Vector2d relative_velocity = point_point_relative_velocity(dp0, dp1);
    CHECK((dp0 - relative_velocity).norm() == Catch::Approx(0).margin(1e-12));
}

TEST_CASE(
    "Point-point relative velocity matrix Jacobians",
    "[friction][point-point][relative_velocity][jacobian]")
{
    const int dim = GENERATE(2, 3);

    auto J_analytical = point_point_relative_velocity_dx_dbeta(dim);

    Eigen::VectorXd x = Eigen::VectorXd::Zero(1);
    Eigen::MatrixXd J_numerical;
    fd::finite_jacobian_tensor<4>(
        x,
        [&dim](const Eigen::VectorXd&) -> Eigen::MatrixXd {
            return point_point_relative_velocity_jacobian(dim);
        },
        J_numerical);

    CHECK(fd::compare_jacobian(J_analytical, J_numerical));
    if (!fd::compare_jacobian(J_analytical, J_numerical)) {
        std::cout << "Analytical:\n"
                  << J_analytical << "\nNumerical:\n"
                  << J_numerical << '\n';
    }
}

TEST_CASE(
    "Point-edge relative velocity matrix Jacobians",
    "[friction][point-edge][relative_velocity][jacobian]")
{
    const int dim = GENERATE(2, 3);
    Eigen::VectorXd x = Eigen::VectorXd::Random(1);

    auto J_analytical = point_edge_relative_velocity_dx_dbeta(dim, x[0]);

    Eigen::MatrixXd J_numerical;
    fd::finite_jacobian_tensor<4>(
        x,
        [&dim](const Eigen::VectorXd& _x) -> Eigen::MatrixXd {
            return point_edge_relative_velocity_jacobian(dim, _x[0]);
        },
        J_numerical);

    CHECK(fd::compare_jacobian(J_analytical, J_numerical));
    if (!fd::compare_jacobian(J_analytical, J_numerical)) {
        std::cout << "Analytical:\n"
                  << J_analytical << "\nNumerical:\n"
                  << J_numerical << '\n';
    }
}

TEST_CASE(
    "Edge-edge relative velocity matrix Jacobians",
    "[friction][edge-edge][relative_velocity][jacobian]")
{
    Eigen::VectorXd x = Eigen::VectorXd::Random(2);

    auto J_analytical = edge_edge_relative_velocity_dx_dbeta(x);

    Eigen::MatrixXd J_numerical;
    fd::finite_jacobian_tensor<4>(
        x,
        [](const Eigen::VectorXd& _x) -> Eigen::MatrixXd {
            return edge_edge_relative_velocity_jacobian(_x);
        },
        J_numerical);

    CHECK(fd::compare_jacobian(J_analytical, J_numerical));
    if (!fd::compare_jacobian(J_analytical, J_numerical)) {
        std::cout << "Analytical:\n"
                  << J_analytical << "\nNumerical:\n"
                  << J_numerical << '\n';
    }
}

TEST_CASE(
    "Point-triangle relative velocity matrix Jacobians",
    "[friction][point-triangle][relative_velocity][jacobian]")
{
    Eigen::VectorXd x = Eigen::VectorXd::Random(2);

    auto J_analytical = point_triangle_relative_velocity_dx_dbeta(x);

    Eigen::MatrixXd J_numerical;
    fd::finite_jacobian_tensor<4>(
        x,
        [](const Eigen::VectorXd& _x) -> Eigen::MatrixXd {
            return point_triangle_relative_velocity_jacobian(_x);
        },
        J_numerical);

    CHECK(fd::compare_jacobian(J_analytical, J_numerical));
    if (!fd::compare_jacobian(J_analytical, J_numerical)) {
        std::cout << "Analytical:\n"
                  << J_analytical << "\nNumerical:\n"
                  << J_numerical << '\n';
    }
}