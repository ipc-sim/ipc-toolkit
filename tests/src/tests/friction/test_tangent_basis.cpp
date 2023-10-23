#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <ipc/friction/tangent_basis.hpp>

#include <finitediff.hpp>

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
    const Eigen::Matrix<double, 36, 2> J =
        point_triangle_tangent_basis_jacobian(p, t0, t1, t2);

    Vector12d x;
    x << p, t0, t1, t2;

    Eigen::MatrixXd tmp;
    fd::finite_jacobian(
        x,
        [](const Eigen::VectorXd& _x) -> Eigen::VectorXd {
            return point_triangle_tangent_basis(
                       _x.segment<3>(0), _x.segment<3>(3), _x.segment<3>(6),
                       _x.segment<3>(9))
                .reshaped();
        },
        tmp);

    Eigen::Matrix<double, 36, 2> J_fd;
    J_fd.col(0) = tmp.topRows(3).reshaped();
    J_fd.col(1) = tmp.bottomRows(3).reshaped();

    CHECK(fd::compare_jacobian(J, J_fd));
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
    const Eigen::Matrix<double, 36, 2> J =
        edge_edge_tangent_basis_jacobian(ea0, ea1, eb0, eb1);

    Vector12d x;
    x << ea0, ea1, eb0, eb1;

    Eigen::MatrixXd tmp;
    fd::finite_jacobian(
        x,
        [](const Eigen::VectorXd& _x) -> Eigen::VectorXd {
            return edge_edge_tangent_basis(
                       _x.segment<3>(0), _x.segment<3>(3), _x.segment<3>(6),
                       _x.segment<3>(9))
                .reshaped();
        },
        tmp);

    Eigen::Matrix<double, 36, 2> J_fd;
    J_fd.col(0) = tmp.topRows(3).reshaped();
    J_fd.col(1) = tmp.bottomRows(3).reshaped();

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
    const Eigen::Matrix<double, 27, 2> J =
        point_edge_tangent_basis_jacobian(p, e0, e1);

    Vector9d x;
    x << p, e0, e1;

    Eigen::MatrixXd tmp;
    fd::finite_jacobian(
        x,
        [](const Eigen::VectorXd& _x) -> Eigen::VectorXd {
            return point_edge_tangent_basis(
                       _x.segment<3>(0), _x.segment<3>(3), _x.segment<3>(6))
                .reshaped();
        },
        tmp);

    Eigen::Matrix<double, 27, 2> J_fd;
    J_fd.col(0) = tmp.topRows(3).reshaped();
    J_fd.col(1) = tmp.bottomRows(3).reshaped();

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
    const Eigen::Matrix<double, 18, 2> J =
        point_point_tangent_basis_jacobian(p0, p1);

    Vector6d x;
    x << p0, p1;

    Eigen::MatrixXd tmp;
    fd::finite_jacobian(
        x,
        [](const Eigen::VectorXd& _x) -> Eigen::VectorXd {
            return point_point_tangent_basis(_x.head<3>(), _x.tail<3>())
                .reshaped();
        },
        tmp);

    Eigen::Matrix<double, 18, 2> J_fd;
    J_fd.col(0) = tmp.topRows(3).reshaped();
    J_fd.col(1) = tmp.bottomRows(3).reshaped();

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
    const Vector12d J = point_edge_tangent_basis_jacobian(p, e0, e1);

    Vector6d x;
    x << p, e0, e1;

    Eigen::MatrixXd J_fd;
    fd::finite_jacobian(
        x,
        [](const Eigen::VectorXd& _x) -> Eigen::VectorXd {
            return point_edge_tangent_basis(
                _x.segment<2>(0), _x.segment<2>(2), _x.segment<2>(4));
        },
        J_fd);
    J_fd = J_fd.reshaped();

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
    const Eigen::Vector<double, 8> J =
        point_point_tangent_basis_jacobian(p0, p1);

    Eigen::Vector4d x;
    x << p0, p1;

    Eigen::MatrixXd J_fd;
    fd::finite_jacobian(
        x,
        [](const Eigen::VectorXd& _x) -> Eigen::VectorXd {
            return point_point_tangent_basis(_x.head<2>(), _x.tail<2>());
        },
        J_fd);
    J_fd = J_fd.reshaped();

    CHECK(fd::compare_jacobian(J, J_fd));
}
