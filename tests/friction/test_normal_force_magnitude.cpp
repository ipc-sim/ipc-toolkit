#include <catch2/catch.hpp>

#include <ipc/friction/normal_force_magnitude.hpp>

#include <ipc/distance/point_triangle.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/point_point.hpp>

#include <finitediff.hpp>

using namespace ipc;

TEST_CASE(
    "Point-triangle normal force magnitude",
    "[friction][point-triangle][normal_force_magnitude]")
{
    Eigen::Vector3d p(0, 1e-4, 0), t0(-1, 0, 1), t1(1, 0, 1), t2(0, 0, -1);

    const double dhat = 1e-3, barrier_stiffness = 1e2;

    const double distance = point_triangle_distance(p, t0, t1, t2);
    VectorMax12d distance_grad;
    point_triangle_distance_gradient(p, t0, t1, t2, distance_grad);
    const VectorMax12d grad = compute_normal_force_magnitude_gradient(
        distance, distance_grad, dhat, barrier_stiffness);

    auto N = [dhat, barrier_stiffness](const Eigen::VectorXd& x) {
        const double distance = point_triangle_distance(
            x.segment<3>(0), x.segment<3>(3), x.segment<3>(6), x.segment<3>(9));
        return compute_normal_force_magnitude(
            distance, dhat, barrier_stiffness);
    };

    Eigen::VectorXd x(12);
    x << p, t0, t1, t2;

    Eigen::VectorXd fgrad;
    fd::finite_gradient(x, N, fgrad);
    CHECK(fd::compare_gradient(grad, fgrad));
}

TEST_CASE(
    "Edge-edge normal force magnitude",
    "[friction][point-triangle][normal_force_magnitude]")
{
    Eigen::Vector3d ea0(-1, -1e-4, 0), ea1(1, -1e-4, 0);
    Eigen::Vector3d eb0(0, 1e-4, -1), eb1(0, 1e-4, 1);

    const double dhat = 1e-3, barrier_stiffness = 1e2;

    const double distance = edge_edge_distance(ea0, ea1, eb0, eb1);
    VectorMax12d distance_grad;
    edge_edge_distance_gradient(ea0, ea1, eb0, eb1, distance_grad);
    const VectorMax12d grad = compute_normal_force_magnitude_gradient(
        distance, distance_grad, dhat, barrier_stiffness);

    auto N = [dhat, barrier_stiffness](const Eigen::VectorXd& x) {
        const double distance = edge_edge_distance(
            x.segment<3>(0), x.segment<3>(3), x.segment<3>(6), x.segment<3>(9));
        return compute_normal_force_magnitude(
            distance, dhat, barrier_stiffness);
    };

    Eigen::VectorXd x(12);
    x << ea0, ea1, eb0, eb1;

    Eigen::VectorXd fgrad;
    fd::finite_gradient(x, N, fgrad);
    CHECK(fd::compare_gradient(grad, fgrad));
}

TEST_CASE(
    "Point-edge normal force magnitude",
    "[friction][point-triangle][normal_force_magnitude]")
{
    Eigen::Vector3d p(0, 1e-4, 0), e0(-1, 0, 0), e1(1, 0, 0);

    const double dhat = 1e-3, barrier_stiffness = 1e2;

    const double distance = point_edge_distance(p, e0, e1);
    VectorMax12d distance_grad;
    point_edge_distance_gradient(p, e0, e1, distance_grad);
    const VectorMax12d grad = compute_normal_force_magnitude_gradient(
        distance, distance_grad, dhat, barrier_stiffness);

    auto N = [dhat, barrier_stiffness](const Eigen::VectorXd& x) {
        const double distance = point_edge_distance(
            x.segment<3>(0), x.segment<3>(3), x.segment<3>(6));
        return compute_normal_force_magnitude(
            distance, dhat, barrier_stiffness);
    };

    Eigen::VectorXd x(9);
    x << p, e0, e1;

    Eigen::VectorXd fgrad;
    fd::finite_gradient(x, N, fgrad);
    CHECK(fd::compare_gradient(grad, fgrad));
}

TEST_CASE(
    "Point-point normal force magnitude",
    "[friction][point-triangle][normal_force_magnitude]")
{
    Eigen::Vector3d p0(0, 0, 0), p1(0, 0, 1e-4);

    const double dhat = 1e-3, barrier_stiffness = 1e2;

    const double distance = point_point_distance(p0, p1);
    VectorMax12d distance_grad;
    point_point_distance_gradient(p0, p1, distance_grad);
    const VectorMax12d grad = compute_normal_force_magnitude_gradient(
        distance, distance_grad, dhat, barrier_stiffness);

    auto N = [dhat, barrier_stiffness](const Eigen::VectorXd& x) {
        const double distance = point_point_distance(x.head<3>(), x.tail<3>());
        return compute_normal_force_magnitude(
            distance, dhat, barrier_stiffness);
    };

    Eigen::VectorXd x(6);
    x << p0, p1;

    Eigen::VectorXd fgrad;
    fd::finite_gradient(x, N, fgrad);
    CHECK(fd::compare_gradient(grad, fgrad));
}

TEST_CASE(
    "Point-edge normal force magnitude (2D)",
    "[friction][point-triangle][normal_force_magnitude]")
{
    Eigen::Vector2d p(0, 1e-4), e0(-1, 0), e1(1, 0);
    const double dhat = 1e-3, barrier_stiffness = 1e2;

    const double distance = point_edge_distance(p, e0, e1);
    VectorMax12d distance_grad;
    point_edge_distance_gradient(p, e0, e1, distance_grad);
    const VectorMax12d grad = compute_normal_force_magnitude_gradient(
        distance, distance_grad, dhat, barrier_stiffness);

    auto N = [dhat, barrier_stiffness](const Eigen::VectorXd& x) {
        const double distance = point_edge_distance(
            x.segment<2>(0), x.segment<2>(2), x.segment<2>(4));
        return compute_normal_force_magnitude(
            distance, dhat, barrier_stiffness);
    };

    Eigen::VectorXd x(6);
    x << p, e0, e1;

    Eigen::VectorXd fgrad;
    fd::finite_gradient(x, N, fgrad);
    CHECK(fd::compare_gradient(grad, fgrad));
}

TEST_CASE(
    "Point-point normal force magnitude (2D)",
    "[friction][point-triangle][normal_force_magnitude][2D]")
{
    Eigen::Vector2d p0(0, 0), p1(0, 1e-4);
    const double dhat = 1e-3, barrier_stiffness = 1e2;

    const double distance = point_point_distance(p0, p1);
    VectorMax12d distance_grad;
    point_point_distance_gradient(p0, p1, distance_grad);
    const VectorMax12d grad = compute_normal_force_magnitude_gradient(
        distance, distance_grad, dhat, barrier_stiffness);

    auto N = [dhat, barrier_stiffness](const Eigen::VectorXd& x) {
        const double distance = point_point_distance(x.head<2>(), x.tail<2>());
        return compute_normal_force_magnitude(
            distance, dhat, barrier_stiffness);
    };

    Eigen::VectorXd x(4);
    x << p0, p1;

    Eigen::VectorXd fgrad;
    fd::finite_gradient(x, N, fgrad);
    CHECK(fd::compare_gradient(grad, fgrad));
}