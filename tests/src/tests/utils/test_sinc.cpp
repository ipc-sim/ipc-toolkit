#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>

#include <iostream>

#include <finitediff.hpp>
#include <igl/PI.h>

// #include <autodiff/autodiff_types.hpp>
#include <ipc/utils/sinc.hpp>

using namespace ipc;

TEST_CASE("sinc(double)", "[sinc][double]")
{
    CHECK(sinc(0) == Catch::Approx(1));
    CHECK(sinc(1e-8) == Catch::Approx(1));
    CHECK(sinc(igl::PI / 2) == Catch::Approx(2 / igl::PI));
    CHECK(sinc(igl::PI) == Catch::Approx(0).margin(1e-16));
}

TEST_CASE("sinc(Interval)", "[sinc][interval]")
{
    // All of these cases fall within the monotonic bound on sinc, so they are
    // a tight bound (ignoring rounding).
    filib::Interval x, expected_y;

    SECTION("y=[1, 1]")
    {
        SECTION("[0, 0]") { x = filib::Interval(0); };
        SECTION("[0, 1e-8]") { x = filib::Interval(0, 1e-8); };
        SECTION("[-1e-8, 0]") { x = filib::Interval(-1e-8, 0); };
        SECTION("[-1e-8, 1e-8]") { x = filib::Interval(-1e-8, 1e-8); };
        expected_y = filib::Interval(1);
    }
    SECTION("y=[0, 1]")
    {
        SECTION("[0, π]") { x = filib::Interval(0, igl::PI); }
        SECTION("[-π, 0]") { x = filib::Interval(-igl::PI, 0); }
        SECTION("[-π, π]") { x = filib::Interval(-igl::PI, igl::PI); }
        expected_y = filib::Interval(0, 1);
    }
    SECTION("y=[2/pi, 1]")
    {
        const double PI_2 = igl::PI / 2;
        SECTION("[0, π/2]") { x = filib::Interval(0, PI_2); }
        SECTION("[-π/2, 0]") { x = filib::Interval(-PI_2, 0); }
        SECTION("[-π/2, π/2]") { x = filib::Interval(-PI_2, PI_2); }
        expected_y = filib::Interval(2 / igl::PI, 1);
    }

    CAPTURE(x.INF, x.SUP);
    filib::Interval y = sinc(x);
    CHECK(y.INF == Catch::Approx(expected_y.INF).margin(1e-8));
    CHECK(y.SUP == Catch::Approx(expected_y.SUP).margin(1e-8));
}

TEST_CASE("Interval sinc with looser bounds", "[sinc][interval]")
{
    filib::Interval x, expected_y;

    SECTION("Non-monotonic")
    {
        const double a = GENERATE(-5, -1, 0);
        x = filib::Interval(a, 5);
        expected_y = filib::Interval(-0.217233628211221659, 1);
    }
    SECTION("Monotonic outside bounds")
    {
        x = filib::Interval(5, 7);
        expected_y = filib::Interval(sinc(5), sinc(7));
    }
    SECTION("Monotonic far outside bounds")
    {
        x = filib::Interval(21, 22);
        expected_y = filib::Interval(sinc(22), sinc(21));
    }

    filib::Interval y = sinc(x);
    CAPTURE(x.INF, x.SUP, y.INF, y.SUP, expected_y.INF, expected_y.SUP);
    CHECK(in(expected_y, y));
    // Loosest bound
    CHECK(in(y, filib::Interval(-0.23, 1)));
    if (x.INF > 0) {
        // Tighter bound if x.INF > 1
        CHECK(in(y, filib::Interval(-1 / x.INF, 1 / x.INF)));
    }
}

TEST_CASE("Interval sinc_norm_x", "[sinc][interval]")
{
    VectorMax3<filib::Interval> x = VectorMax3<filib::Interval>::Zero(3);
    filib::Interval expected_y;

    SECTION("Zero")
    {
        // x is already zero
        expected_y = filib::Interval(1);
    }
    SECTION("Positive")
    {
        x(1) = filib::Interval(igl::PI);
        expected_y = filib::Interval(0);
    }
    SECTION("Negative")
    {
        x(1) = filib::Interval(-igl::PI);
        expected_y = filib::Interval(0);
    }
    SECTION("Mixed")
    {
        x(1) = filib::Interval(-igl::PI, igl::PI);
        expected_y = filib::Interval(0, 1);
    }

    filib::Interval y = sinc_norm_x<filib::Interval>(x);
    CAPTURE(y.INF, y.SUP, expected_y.INF, expected_y.SUP);
    CHECK(expected_y.INF == Catch::Approx(y.INF).margin(1e-8));
    CHECK(expected_y.SUP == Catch::Approx(y.SUP).margin(1e-8));
}

TEST_CASE("∇sinc(||x||)", "[sinc][vector][diff]")
{
    double sign = GENERATE(-1, 1);
    double val = GENERATE(0, 1e-8, igl::PI / 2, igl::PI, 5, 20, 100);
    int index = GENERATE(0, 1, 2);
    Eigen::Vector3d x = Eigen::Vector3d::Zero();
    x(index) = sign * val;

    Eigen::VectorXd fgrad(3);
    fd::finite_gradient(x, sinc_norm_x<double>, fgrad);

    Eigen::Vector3d grad = sinc_norm_x_grad(x);
    CHECK(fd::compare_gradient(grad, fgrad));
}

TEST_CASE("∇²sinc(||x||)", "[sinc][vector][diff]")
{
    double sign = GENERATE(-1, 1);
    double val = GENERATE(0, 1e-8, igl::PI / 2, igl::PI, 5, 20, 100);
    int index = GENERATE(0, 1, 2);
    Eigen::Vector3d x = Eigen::Vector3d::Zero();
    x(index) = sign * val;

    Eigen::MatrixXd fhess(3, 3);
    fd::finite_hessian(x, sinc_norm_x<double>, fhess);

    Eigen::Matrix3d hess = sinc_norm_x_hess(x);
    CHECK(fd::compare_hessian(hess, fhess));
}