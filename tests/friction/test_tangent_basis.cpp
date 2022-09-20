#include <catch2/catch.hpp>

#include <ipc/friction/tangent_basis.hpp>

using namespace ipc;

TEST_CASE(
    "Point-triangle tangent basis", "[friction][point-triangle][tangent_basis]")
{
    Eigen::Vector3d p(0, 1, 0), t0(-1, 0, 1), t1(1, 0, 1), t2(0, 0, -1);

    Eigen::Matrix<double, 3, 2> basis =
        point_triangle_tangent_basis(p, t0, t1, t2);
    CAPTURE(basis);
    CHECK(std::abs(basis.col(0).dot(Eigen::Vector3d::UnitX())) == Approx(1));
    CHECK(std::abs(basis.col(1).dot(Eigen::Vector3d::UnitZ())) == Approx(1));
}

TEST_CASE("Edge-edge tangent basis", "[friction][edge-edge][tangent_basis]")
{
    Eigen::Vector3d ea0(-1, 0, 0), ea1(1, 0, 0), eb0(0, 0, -1), eb1(0, 0, 1);

    Eigen::Matrix<double, 3, 2> basis =
        edge_edge_tangent_basis(ea0, ea1, eb0, eb1);
    CAPTURE(basis);
    CHECK(std::abs(basis.col(0).dot(Eigen::Vector3d::UnitX())) == Approx(1));
    CHECK(std::abs(basis.col(1).dot(Eigen::Vector3d::UnitZ())) == Approx(1));
}

TEST_CASE("Point-edge tangent basis", "[friction][point-edge][tangent_basis]")
{
    Eigen::Vector3d p(0, 1, 0), e0(-1, 0, 0), e1(1, 0, 0);

    Eigen::Matrix<double, 3, 2> basis = point_edge_tangent_basis(p, e0, e1);
    CAPTURE(basis);
    CHECK(std::abs(basis.col(0).dot(Eigen::Vector3d::UnitX())) == Approx(1));
    CHECK(std::abs(basis.col(1).dot(Eigen::Vector3d::UnitZ())) == Approx(1));
}

TEST_CASE("Point-point tangent basis", "[friction][point-point][tangent_basis]")
{
    Eigen::Vector3d p0(0, 0, 0), p1(0, 0, 1);

    Eigen::Matrix<double, 3, 2> basis = point_point_tangent_basis(p0, p1);
    CAPTURE(basis);
    CHECK(std::abs(basis.col(0).dot(Eigen::Vector3d::UnitX())) == Approx(1));
    CHECK(std::abs(basis.col(1).dot(Eigen::Vector3d::UnitY())) == Approx(1));
}

TEST_CASE(
    "Point-edge tangent basis in 2D",
    "[friction][point-edge][tangent_basis][2D]")
{
    Eigen::Vector2d p(0, 1), e0(-1, 0), e1(1, 0);

    Eigen::Matrix<double, 2, 1> basis = point_edge_tangent_basis(p, e0, e1);
    CAPTURE(basis);
    CHECK(std::abs(basis.dot(Eigen::Vector2d::UnitX())) == Approx(1));
}

TEST_CASE(
    "Point-point tangent basis in 2D",
    "[friction][point-point][tangent_basis][2D]")
{
    Eigen::Vector2d p0(0, 0), p1(0, 1);

    Eigen::Matrix<double, 2, 1> basis = point_point_tangent_basis(p0, p1);
    CAPTURE(basis);
    CHECK(std::abs(basis.dot(Eigen::Vector2d::UnitX())) == Approx(1));
}
