#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <ipc/tangent/tangent_basis.hpp>

#include <finitediff.hpp>

#include <iostream>

using namespace ipc;

TEST_CASE(
    "Point-triangle tangent basis",
    "[friction][point-triangle][tangent_basis][tangent_basis_jacobian]")
{
    Eigen::Vector3d p(0, 1, 0), t0(-1, 0, 1), t1(1, 0, 1), t2(0, 0, -1);

    Eigen::Matrix<double, 3, 2> basis =
        point_triangle_tangent_basis(p, t0, t1, t2);
    CAPTURE(basis);
    CHECK(
        std::abs(basis.col(0).dot(Eigen::Vector3d::UnitX()))
        == Catch::Approx(1));
    CHECK(
        std::abs(basis.col(1).dot(Eigen::Vector3d::UnitZ()))
        == Catch::Approx(1));

    // Jacobian
    const auto J = point_triangle_tangent_basis_jacobian(p, t0, t1, t2);

    Vector12d x;
    x << p, t0, t1, t2;

    Eigen::MatrixXd J_fd;
    fd::finite_jacobian_tensor<3>(
        x,
        [](const Eigen::VectorXd& _x) -> Eigen::MatrixXd {
            return point_triangle_tangent_basis(
                _x.segment<3>(0), _x.segment<3>(3), _x.segment<3>(6),
                _x.segment<3>(9));
        },
        J_fd);

    CHECK(fd::compare_jacobian(J, J_fd));
    if (!fd::compare_jacobian(J, J_fd)) {
        std::cout << "J_analytic:\n" << J << "\n\n";
        std::cout << "J_numerical:\n" << J_fd << "\n\n";
    }
}

TEST_CASE(
    "Edge-edge tangent basis",
    "[friction][edge-edge][tangent_basis][tangent_basis_jacobian]")
{
    Eigen::Vector3d ea0(-1, 0, 0), ea1(1, 0, 0), eb0(0, 0, -1), eb1(0, 0, 1);

    Eigen::Matrix<double, 3, 2> basis =
        edge_edge_tangent_basis(ea0, ea1, eb0, eb1);
    CAPTURE(basis);
    CHECK(
        std::abs(basis.col(0).dot(Eigen::Vector3d::UnitX()))
        == Catch::Approx(1));
    CHECK(
        std::abs(basis.col(1).dot(Eigen::Vector3d::UnitZ()))
        == Catch::Approx(1));

    // Jacobian
    const auto J = edge_edge_tangent_basis_jacobian(ea0, ea1, eb0, eb1);

    Vector12d x;
    x << ea0, ea1, eb0, eb1;

    Eigen::MatrixXd J_fd;
    fd::finite_jacobian_tensor<3>(
        x,
        [](const Eigen::VectorXd& _x) -> Eigen::MatrixXd {
            return edge_edge_tangent_basis(
                _x.segment<3>(0), _x.segment<3>(3), _x.segment<3>(6),
                _x.segment<3>(9));
        },
        J_fd);

    CHECK(fd::compare_jacobian(J, J_fd));
}

TEST_CASE(
    "Point-edge tangent basis",
    "[friction][point-edge][tangent_basis][tangent_basis_jacobian]")
{
    Eigen::Vector3d p(0, 1, 0), e0(-1, 0, 0), e1(1, 0, 0);

    Eigen::Matrix<double, 3, 2> basis = point_edge_tangent_basis(p, e0, e1);
    CAPTURE(basis);
    CHECK(
        std::abs(basis.col(0).dot(Eigen::Vector3d::UnitX()))
        == Catch::Approx(1));
    CHECK(
        std::abs(basis.col(1).dot(Eigen::Vector3d::UnitZ()))
        == Catch::Approx(1));

    // Jacobian
    const auto J = point_edge_tangent_basis_jacobian(p, e0, e1);

    Vector9d x;
    x << p, e0, e1;

    Eigen::MatrixXd J_fd;
    fd::finite_jacobian_tensor<3>(
        x,
        [](const Eigen::VectorXd& _x) -> Eigen::MatrixXd {
            return point_edge_tangent_basis(
                _x.segment<3>(0), _x.segment<3>(3), _x.segment<3>(6));
        },
        J_fd);

    CHECK(fd::compare_jacobian(J, J_fd));
}

TEST_CASE(
    "Point-point tangent basis",
    "[friction][point-point][tangent_basis][tangent_basis_jacobian]")
{
    Eigen::Vector3d p0(0, 0, 0), p1(0, 0, 1);

    Eigen::Matrix<double, 3, 2> basis = point_point_tangent_basis(p0, p1);
    CAPTURE(basis);
    CHECK(
        std::abs(basis.col(0).dot(Eigen::Vector3d::UnitX()))
        == Catch::Approx(1));
    CHECK(
        std::abs(basis.col(1).dot(Eigen::Vector3d::UnitY()))
        == Catch::Approx(1));

    // Jacobian
    const auto J = point_point_tangent_basis_jacobian(p0, p1);

    Vector6d x;
    x << p0, p1;

    Eigen::MatrixXd J_fd;
    fd::finite_jacobian_tensor<3>(
        x,
        [](const Eigen::VectorXd& _x) -> Eigen::MatrixXd {
            return point_point_tangent_basis(_x.head<3>(), _x.tail<3>());
        },
        J_fd);

    CHECK(fd::compare_jacobian(J, J_fd));
}

TEST_CASE(
    "Point-edge tangent basis in 2D",
    "[friction][point-edge][tangent_basis][tangent_basis_jacobian][2D]")
{
    Eigen::Vector2d p(0, 1), e0(-1, 0), e1(1, 0);

    Eigen::Matrix<double, 2, 1> basis = point_edge_tangent_basis(p, e0, e1);
    CAPTURE(basis);
    CHECK(std::abs(basis.dot(Eigen::Vector2d::UnitX())) == Catch::Approx(1));

    // Jacobian
    const auto J = point_edge_tangent_basis_jacobian(p, e0, e1);

    Vector6d x;
    x << p, e0, e1;

    Eigen::MatrixXd J_fd;
    fd::finite_jacobian_tensor<3>(
        x,
        [](const Eigen::VectorXd& _x) -> Eigen::VectorXd {
            return point_edge_tangent_basis(
                _x.segment<2>(0), _x.segment<2>(2), _x.segment<2>(4));
        },
        J_fd);

    CHECK(fd::compare_jacobian(J, J_fd));
}

TEST_CASE(
    "Point-point tangent basis in 2D",
    "[friction][point-point][tangent_basis][tangent_basis_jacobian][2D]")
{
    Eigen::Vector2d p0(0, 0), p1(0, 1);

    Eigen::Matrix<double, 2, 1> basis = point_point_tangent_basis(p0, p1);
    CAPTURE(basis);
    CHECK(std::abs(basis.dot(Eigen::Vector2d::UnitX())) == Catch::Approx(1));

    // Jacobian
    const auto J = point_point_tangent_basis_jacobian(p0, p1);

    Eigen::Vector4d x;
    x << p0, p1;

    Eigen::MatrixXd J_fd;
    fd::finite_jacobian_tensor<3>(
        x,
        [](const Eigen::VectorXd& _x) -> Eigen::VectorXd {
            return point_point_tangent_basis(_x.head<2>(), _x.tail<2>());
        },
        J_fd);

    CHECK(fd::compare_jacobian(J, J_fd));
}

TEST_CASE(
    "Tangent bases agree whether or not the dimension is known at compile time",
    "[friction][tangent_basis][tangent_basis_jacobian]")
{
    // The front ends dispatch on `dim_v<Derived>` when the argument type
    // carries its size, and fall back to a runtime branch on `size()` when it
    // does not. Every other test here passes a fixed-size vector, so only the
    // `if constexpr` arms ever run. Passing the same points as a dynamically
    // sized vector takes the fallback instead, and the two must agree exactly:
    // they call the same kernel, so any difference is a dispatch bug.
    const int dim = GENERATE(2, 3);
    CAPTURE(dim);

    const Eigen::VectorXd p = Eigen::VectorXd::LinSpaced(dim, 0.25, 1.0);
    const Eigen::VectorXd q = Eigen::VectorXd::LinSpaced(dim, -1.0, 0.5);

    // VectorMax3d keeps its size at runtime, so dim_v is not a constant.
    const VectorMax3d p_dyn = p, q_dyn = q;

    const Eigen::VectorXd r = Eigen::VectorXd::LinSpaced(dim, 0.75, -0.5);
    const VectorMax3d r_dyn = r;

    if (dim == 2) {
        const Eigen::Vector2d a = p, b = q, c = r;
        CHECK(
            point_point_tangent_basis(p_dyn, q_dyn)
            == point_point_tangent_basis(a, b));
        CHECK(
            point_point_tangent_basis_jacobian(p_dyn, q_dyn)
            == point_point_tangent_basis_jacobian(a, b));
        CHECK(
            point_edge_tangent_basis(p_dyn, q_dyn, r_dyn)
            == point_edge_tangent_basis(a, b, c));
        CHECK(
            point_edge_tangent_basis_jacobian(p_dyn, q_dyn, r_dyn)
            == point_edge_tangent_basis_jacobian(a, b, c));
    } else {
        const Eigen::Vector3d a = p, b = q, c = r;
        CHECK(
            point_point_tangent_basis(p_dyn, q_dyn)
            == point_point_tangent_basis(a, b));
        CHECK(
            point_point_tangent_basis_jacobian(p_dyn, q_dyn)
            == point_point_tangent_basis_jacobian(a, b));
        CHECK(
            point_edge_tangent_basis(p_dyn, q_dyn, r_dyn)
            == point_edge_tangent_basis(a, b, c));
        CHECK(
            point_edge_tangent_basis_jacobian(p_dyn, q_dyn, r_dyn)
            == point_edge_tangent_basis_jacobian(a, b, c));
    }
}
