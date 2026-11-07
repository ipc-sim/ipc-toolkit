#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <ipc/tangent/closest_point.hpp>

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

    // test Hessian: one 12x12 block per barycentric coordinate, each the
    // Jacobian of that coordinate's row of J.
    const std::array<Eigen::Matrix<double, 12, 12>, 2> H =
        point_triangle_closest_point_hessian(p, t0, t1, t2);

    for (int c = 0; c < 2; c++) {
        CAPTURE(c);
        Eigen::MatrixXd H_FD;
        fd::finite_jacobian(
            x,
            [c](const Eigen::VectorXd& _x) -> Eigen::VectorXd {
                return point_triangle_closest_point_jacobian(
                           _x.segment<3>(0), _x.segment<3>(3), _x.segment<3>(6),
                           _x.segment<3>(9))
                    .row(c)
                    .transpose();
            },
            H_FD);

        CHECK(fd::compare_jacobian(H[c], H_FD));
    }
}

TEST_CASE("Edge-edge closest point", "[friction][edge-edge][closest_point]")
{
    Eigen::Vector3d ea0(-1, 0, 0), ea1(1, 0, 0), eb0(0, 0, -1), eb1(0, 0, 1);

    Eigen::Vector2d barycentric_coords =
        edge_edge_closest_point(ea0, ea1, eb0, eb1);
    CAPTURE(barycentric_coords);
    // Perpendicular edges centered on the same axis meet at their midpoints,
    // so the answer is exactly (0.5, 0.5) with no rounding to hide behind.
    CHECK(barycentric_coords[0] == Catch::Approx(0.5).epsilon(0).margin(1e-15));
    CHECK(barycentric_coords[1] == Catch::Approx(0.5).epsilon(0).margin(1e-15));

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

    // test Hessian: one 12x12 block per barycentric coordinate, each the
    // Jacobian of that coordinate's row of J.
    const std::array<Eigen::Matrix<double, 12, 12>, 2> H =
        edge_edge_closest_point_hessian(ea0, ea1, eb0, eb1);

    for (int c = 0; c < 2; c++) {
        CAPTURE(c);
        Eigen::MatrixXd H_FD;
        fd::finite_jacobian(
            x,
            [c](const Eigen::VectorXd& _x) -> Eigen::VectorXd {
                return edge_edge_closest_point_jacobian(
                           _x.segment<3>(0), _x.segment<3>(3), _x.segment<3>(6),
                           _x.segment<3>(9))
                    .row(c)
                    .transpose();
            },
            H_FD);

        CHECK(fd::compare_jacobian(H[c], H_FD));
    }
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

    // test Hessian: the closest point is a scalar here, so its Hessian is the
    // Jacobian of the gradient checked above.
    const Eigen::Matrix<double, 9, 9> H =
        point_edge_closest_point_hessian(p, e0, e1);

    Eigen::MatrixXd H_FD;
    fd::finite_jacobian(
        x,
        [](const Eigen::VectorXd& _x) -> Eigen::VectorXd {
            return point_edge_closest_point_jacobian(
                _x.segment<3>(0), _x.segment<3>(3), _x.segment<3>(6));
        },
        H_FD);

    CHECK(fd::compare_jacobian(H, H_FD));
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

    // test Hessian: covers the 2D branch of the kernel, which the 3D case
    // above never reaches.
    const Eigen::Matrix<double, 6, 6> H =
        point_edge_closest_point_hessian(p, e0, e1);

    Eigen::MatrixXd H_FD;
    fd::finite_jacobian(
        x,
        [](const Eigen::VectorXd& _x) -> Eigen::VectorXd {
            return point_edge_closest_point_jacobian(
                _x.segment<2>(0), _x.segment<2>(2), _x.segment<2>(4));
        },
        H_FD);

    CHECK(fd::compare_jacobian(H, H_FD));
}

TEST_CASE(
    "Point-edge closest point agrees whether or not the dimension is known at "
    "compile time",
    "[friction][point-edge][closest_point]")
{
    // The front ends dispatch on `dim_v<Derived>` when the argument type
    // carries its size and fall back to a runtime branch on `size()` when it
    // does not. Every other test here passes a fixed-size vector, so only the
    // `if constexpr` arms ever run. A VectorMax3d holds the same numbers but
    // keeps its size at runtime, taking the fallback instead; both arms reach
    // the same kernel, so any difference is a dispatch bug.
    const int dim = GENERATE(2, 3);
    CAPTURE(dim);

    const Eigen::VectorXd p = Eigen::VectorXd::LinSpaced(dim, 0.25, 1.0);
    const Eigen::VectorXd e0 = Eigen::VectorXd::LinSpaced(dim, -1.0, 0.5);
    const Eigen::VectorXd e1 = Eigen::VectorXd::LinSpaced(dim, 0.75, -0.5);

    const VectorMax3d p_dyn = p, e0_dyn = e0, e1_dyn = e1;

    if (dim == 2) {
        const Eigen::Vector2d a = p, b = e0, c = e1;
        CHECK(
            point_edge_closest_point(p_dyn, e0_dyn, e1_dyn)
            == point_edge_closest_point(a, b, c));
        CHECK(
            point_edge_closest_point_jacobian(p_dyn, e0_dyn, e1_dyn)
            == point_edge_closest_point_jacobian(a, b, c));
        CHECK(
            point_edge_closest_point_hessian(p_dyn, e0_dyn, e1_dyn)
            == point_edge_closest_point_hessian(a, b, c));
    } else {
        const Eigen::Vector3d a = p, b = e0, c = e1;
        CHECK(
            point_edge_closest_point(p_dyn, e0_dyn, e1_dyn)
            == point_edge_closest_point(a, b, c));
        CHECK(
            point_edge_closest_point_jacobian(p_dyn, e0_dyn, e1_dyn)
            == point_edge_closest_point_jacobian(a, b, c));
        CHECK(
            point_edge_closest_point_hessian(p_dyn, e0_dyn, e1_dyn)
            == point_edge_closest_point_hessian(a, b, c));
    }
}
