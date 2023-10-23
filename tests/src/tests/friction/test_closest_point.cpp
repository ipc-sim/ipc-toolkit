#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <ipc/friction/closest_point.hpp>

#include <finitediff.hpp>

using namespace ipc;

TEST_CASE(
    "Point-triangle closest point", "[friction][point-triangle][closest_point]")
{
    Eigen::Vector3d t0(-1, 0, 1), t1(1, 0, 1), t2(0, 0, -1);
    Eigen::Vector2d expected_coords(0.5, 0.5);
    Eigen::Vector3d p =
        t0 + expected_coords[0] * (t1 - t0) + expected_coords[1] * (t2 - t0);
    // p = 1 * t0 + u * t1 - u * t0 + v * t2 - v * t0
    //   = (1 - u - v) * t0 + u * t1 + v * t2
    //   =  w * t0 + u * t1 + v * t2

    Eigen::Vector2d barycentric_coords =
        point_triangle_closest_point(p, t0, t1, t2);
    Eigen::Vector3d p_actual = t0 + barycentric_coords[0] * (t1 - t0)
        + barycentric_coords[1] * (t2 - t0);
    CAPTURE(barycentric_coords);
    CHECK((p - p_actual).norm() == Catch::Approx(0).margin(1e-12));

    // test Jacobian
    Eigen::Matrix<double, 2, 12> J =
        point_triangle_closest_point_jacobian(p, t0, t1, t2);

    Vector12d x;
    x << p, t0, t1, t2;

    Eigen::MatrixXd J_FD;
    fd::finite_jacobian(
        x,
        [](const Eigen::VectorXd& _x) {
            return point_triangle_closest_point(
                _x.segment<3>(0), _x.segment<3>(3), _x.segment<3>(6),
                _x.segment<3>(9));
        },
        J_FD);

    CHECK(fd::compare_jacobian(J, J_FD));
}

TEST_CASE("Edge-edge closest point", "[friction][edge-edge][closest_point]")
{
    Eigen::Vector3d ea0(-1, 0, 0), ea1(1, 0, 0), eb0(0, 0, -1), eb1(0, 0, 1);

    Eigen::Vector2d barycentric_coords =
        edge_edge_closest_point(ea0, ea1, eb0, eb1);
    CAPTURE(barycentric_coords);
    CHECK(barycentric_coords[0] == Catch::Approx(0.5));
    CHECK(barycentric_coords[1] == Catch::Approx(0.5));

    // test Jacobian
    Eigen::Matrix<double, 2, 12> J =
        edge_edge_closest_point_jacobian(ea0, ea1, eb0, eb1);

    Vector12d x;
    x << ea0, ea1, eb0, eb1;

    Eigen::MatrixXd J_FD;
    fd::finite_jacobian(
        x,
        [](const Eigen::VectorXd& _x) {
            return edge_edge_closest_point(
                _x.segment<3>(0), _x.segment<3>(3), _x.segment<3>(6),
                _x.segment<3>(9));
        },
        J_FD);

    CHECK(fd::compare_jacobian(J, J_FD));
}

TEST_CASE("Point-edge closest point", "[friction][point-edge][closest_point]")
{
    Eigen::Vector3d p(0, 1, 0), e0(-1, 0, 0), e1(1, 0, 0);

    double alpha = point_edge_closest_point(p, e0, e1);
    CHECK(alpha == Catch::Approx(0.5));

    // test Jacobian
    VectorMax9d J = point_edge_closest_point_jacobian(p, e0, e1);

    Vector9d x;
    x << p, e0, e1;

    Eigen::VectorXd J_FD;
    fd::finite_gradient(
        x,
        [](const Eigen::VectorXd& _x) {
            return point_edge_closest_point(
                _x.segment<3>(0), _x.segment<3>(3), _x.segment<3>(6));
        },
        J_FD);

    CHECK(fd::compare_gradient(J, J_FD));
}

TEST_CASE(
    "Point-edge closest point in 2D",
    "[friction][point-edge][closest_point][2D]")
{
    Eigen::Vector2d p(0, 1), e0(-1, 0), e1(1, 0);

    double alpha = point_edge_closest_point(p, e0, e1);
    CHECK(alpha == Catch::Approx(0.5));

    // test Jacobian
    VectorMax9d J = point_edge_closest_point_jacobian(p, e0, e1);

    Vector6d x;
    x << p, e0, e1;

    Eigen::VectorXd J_FD;
    fd::finite_gradient(
        x,
        [](const Eigen::VectorXd& _x) {
            return point_edge_closest_point(
                _x.segment<2>(0), _x.segment<2>(2), _x.segment<2>(4));
        },
        J_FD);

    CHECK(fd::compare_gradient(J, J_FD));
}
