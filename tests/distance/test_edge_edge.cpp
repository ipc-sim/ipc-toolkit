#include <catch2/catch.hpp>

#include <finitediff.hpp>
#include <igl/PI.h>

#include <distance/edge_edge.hpp>
#include <utils/eigen_ext.hpp>

using namespace ipc;

TEST_CASE("Edge-edge distance", "[distance][edge-edge]")
{
    double e0y = GENERATE(-10, -1, -1e-4, 0, 1e-4, 1, 10);
    Eigen::Vector3d e00(-1, e0y, 0);
    Eigen::Vector3d e01(1, e0y, 0);
    Eigen::Vector3d e10(0, 0, -1);
    Eigen::Vector3d e11(0, 0, 1);

    Eigen::Vector3d e0_closest, e1_closest;
    double shiftx = GENERATE(-2, 0, 2);
    double shiftz = GENERATE(-2, 0, 2);
    double e0x = shiftx + GENERATE(-1, -0.5, 0, 0.5, 1);
    double e0z = shiftz + GENERATE(-1, -0.5, 0, 0.5, 1);
    e00.x() += e0x;
    e01.x() += e0x;
    e00.z() += e0z;
    e01.z() += e0z;
    e0_closest =
        shiftx > 1 ? e00 : (shiftx < -1 ? e01 : Eigen::Vector3d(0, e0y, e0z));
    e1_closest =
        shiftz > 1 ? e11 : (shiftz < -1 ? e10 : Eigen::Vector3d(0, 0, e0z));

    double distance = edge_edge_distance(e00, e01, e10, e11);
    double expected_distance = point_point_distance(e0_closest, e1_closest);

    CAPTURE(e0x, e0y, e0z, edge_edge_distance_type(e00, e01, e10, e11));
    CHECK(distance == Approx(expected_distance).margin(1e-12));
}

TEST_CASE("Edge-edge distance degenerate case", "[distance][edge-edge]")
{
    double e0y = GENERATE(-10, -1, -1e-4, 0, 1e-4, 1, 10);
    Eigen::Vector3d e00(-1, e0y, 0);
    Eigen::Vector3d e01(1, e0y, 0);
    Eigen::Vector3d e10(0, 0, -1);
    Eigen::Vector3d e11(0, 0, 1);

    double theta =
        GENERATE(-2, -1.5, -1, -0.123124, 0, 0.2342352, 0.5, 1, 1.5, 2, 50, 51)
        * igl::PI;
    Eigen::Matrix3d R =
        Eigen::AngleAxisd(theta, Eigen::Vector3d::UnitY()).toRotationMatrix();

    SECTION("e0 rotating")
    {
        e00 = R * e00;
        e01 = R * e01;
    }
    SECTION("e1 rotating")
    {
        e10 = R * e10;
        e11 = R * e11;
    }

    double distance = edge_edge_distance(e00, e01, e10, e11);
    CHECK(distance == Approx(e0y * e0y).margin(1e-12));
}

TEST_CASE(
    "Edge-edge distance degenerate case not overlapping",
    "[distance][edge-edge]")
{
    double e0y = GENERATE(-10, -1, -1e-4, 0, 1e-4, 1, 10);
    double gap = GENERATE(0, 0.01, 0.1, 1);
    Eigen::Vector3d e00(gap, e0y, 0);
    Eigen::Vector3d e01(1, e0y, 0);
    Eigen::Vector3d e10(-1, 0, 0);
    Eigen::Vector3d e11(-gap, 0, 0);

    SECTION("original order") {}
    SECTION("swap e0") { std::swap(e00, e01); }
    SECTION("swap e1") { std::swap(e10, e11); }
    SECTION("swap e0 and e1")
    {
        std::swap(e00, e01);
        std::swap(e10, e11);
    }

    double distance = edge_edge_distance(e00, e01, e10, e11);
    double expected_distance = point_point_distance(
        Eigen::Vector3d(gap, e0y, 0), Eigen::Vector3d(-gap, 0, 0));

    CHECK(distance == Approx(expected_distance).margin(1e-12));
}

TEST_CASE("Edge-edge distance gradient", "[distance][edge-edge][gradient]")
{
    double e0y = GENERATE(-10, -1, -1e-4, 0, 1e-4, 1, 10);
    Eigen::Vector3d e00(-1, e0y, 0);
    Eigen::Vector3d e01(1, e0y, 0);
    Eigen::Vector3d e10(0, 0, -1);
    Eigen::Vector3d e11(0, 0, 1);

    Eigen::Vector3d e0_closest, e1_closest;
    double shiftx = GENERATE(-2, 0, 2);
    double shiftz = GENERATE(-2, 0, 2);
    double e0x = shiftx + GENERATE(-1, -0.5, 0, 0.5, 1);
    double e0z = shiftz + GENERATE(-1, -0.5, 0, 0.5, 1);
    e00.x() += e0x;
    e01.x() += e0x;
    e00.z() += e0z;
    e01.z() += e0z;
    e0_closest =
        shiftx > 1 ? e00 : (shiftx < -1 ? e01 : Eigen::Vector3d(0, e0y, e0z));
    e1_closest =
        shiftz > 1 ? e11 : (shiftz < -1 ? e10 : Eigen::Vector3d(0, 0, e0z));

    Eigen::Matrix<double, 12, 1> grad;
    edge_edge_distance_gradient(e00, e01, e10, e11, grad);

    // Compute the gradient using finite differences
    Eigen::VectorXd x(12);
    x.segment<3>(0) = e00;
    x.segment<3>(3) = e01;
    x.segment<3>(6) = e10;
    x.segment<3>(9) = e11;
    auto f = [](const Eigen::VectorXd& x) {
        return edge_edge_distance(
            x.segment<3>(0), x.segment<3>(3), x.segment<3>(6), x.segment<3>(9));
    };
    Eigen::VectorXd fgrad;
    fd::finite_gradient(x, f, fgrad);

    CAPTURE(e0x, e0y, e0z, edge_edge_distance_type(e00, e01, e10, e11));
    CHECK(fd::compare_gradient(grad, fgrad));
}

TEST_CASE(
    "Parallel edge-edge distance gradient", "[distance][edge-edge][gradient]")
{
    // Generate a geometric space of
    double angle = 0;
    SECTION("Almost parallel")
    {
        double exponent = GENERATE(range(-6, 3));
        angle = pow(10, exponent) * igl::PI / 180.0;
    }
    // SECTION("Parallel") { angle = 0; }

    Eigen::Vector3d e00(-1.0, 0, 0), e01(1.0, 0, 0),
        e10(cos(angle), 1, sin(angle)),
        e11(cos(angle + igl::PI), 1, sin(angle + igl::PI));

    double distance = edge_edge_distance(e00, e01, e10, e11);
    Eigen::VectorXd grad;
    edge_edge_distance_gradient(e00, e01, e10, e11, grad);

    // Compute the gradient using finite differences
    Eigen::VectorXd x(12);
    x.segment<3>(0) = e00;
    x.segment<3>(3) = e01;
    x.segment<3>(6) = e10;
    x.segment<3>(9) = e11;
    auto f = [](const Eigen::VectorXd& x) {
        return edge_edge_distance(
            x.segment<3>(0), x.segment<3>(3), x.segment<3>(6), x.segment<3>(9));
    };
    Eigen::VectorXd fgrad;
    fd::finite_gradient(x, f, fgrad);

    CAPTURE(angle, (grad - fgrad).squaredNorm());
    CHECK(distance == Approx(1.0));
    CHECK(fd::compare_gradient(grad, fgrad));
    // CHECK(distance.getHessian().squaredNorm() != Approx(0.0));
}
