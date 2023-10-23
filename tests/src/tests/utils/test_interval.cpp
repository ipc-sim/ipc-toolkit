#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <ipc/config.hpp>

#ifdef IPC_TOOLKIT_WITH_FILIB

#include <ipc/utils/interval.hpp>
#include <ipc/utils/logger.hpp>

#include <igl/PI.h>

using namespace ipc;

TEST_CASE("Simple interval arithmetic", "[interval]")
{
    filib::Interval i(0, 1), j(4, 5);
    CHECK(i.INF <= i.SUP);
    CHECK(j.INF <= j.SUP);

    filib::Interval result;
    result = i + j;
    CHECK(!empty(result));
    result = i - j;
    CHECK(!empty(result));
    result = i * j;
    CHECK(!empty(result));
    result = 1.0 / j;
    CHECK(!empty(result));
    result = i / j;
    CHECK(!empty(result));

    result = i + 10.0;
    CHECK(result.INF == Catch::Approx(i.INF + 10).margin(1e-12));
    CHECK(result.SUP == Catch::Approx(i.SUP + 10).margin(1e-12));
}

TEST_CASE("Cosine interval arithmetic", "[interval]")
{
    filib::Interval r;

    double shift;
    SECTION("No shift") { shift = 0; }
    SECTION("2π shift") { shift = 2 * igl::PI; }
    SECTION("-2π shift") { shift = -2 * igl::PI; }
    SECTION("100π shift") { shift = 100 * igl::PI; }
    SECTION("-100π shift") { shift = -100 * igl::PI; }

    CAPTURE(shift);

    r = cos(filib::Interval(-1, 7) + shift);
    CHECK(r.INF == -1.0);
    CHECK(r.SUP == 1.0);

    r = cos(filib::Interval(2, 4) + shift);
    CHECK(r.SUP < 0);

    r = cos(filib::Interval(0, 1) + shift);
    CHECK(r.INF > 0);

    r = cos(filib::Interval(1, 2) + shift);
    CHECK(r.INF < 0);
    CHECK(r.SUP > 0);
}

TEST_CASE("Sine interval arithmetic", "[interval]")
{
    filib::Interval r;

    double shift = igl::PI / 2;
    SECTION("No shift") { shift += 0; }
    SECTION("2π shift") { shift += 2 * igl::PI; }
    SECTION("-2π shift") { shift += -2 * igl::PI; }
    SECTION("100π shift") { shift += 100 * igl::PI; }
    SECTION("-100π shift") { shift += -100 * igl::PI; }

    CAPTURE(shift);

    r = sin(filib::Interval(-1, 7) + shift);
    CHECK(r.INF == -1.0);
    CHECK(r.SUP == 1.0);

    r = sin(filib::Interval(2, 4) + shift);
    CAPTURE(
        (filib::Interval(2, 4) + shift).INF,
        (filib::Interval(2, 4) + shift).SUP);
    CHECK(r.SUP < 0);

    r = sin(filib::Interval(0, 1) + shift);
    CHECK(r.INF > 0);

    r = sin(filib::Interval(1, 2) + shift);
    CHECK(r.INF < 0);
    CHECK(r.SUP > 0);
}

TEST_CASE("Interval rotation rounding", "[interval][matrix]")
{
    const filib::Interval x(0.30969396267858817);
    const filib::Interval y(0.85675409103416755);
    const filib::Interval theta(0.79358805865013693);

    Eigen::Matrix<double, 4, 2> V;
    V.row(0) << 0.5, 0.0;
    V.row(1) << 0.0, 0.5;
    V.row(2) << -0.5, 0.0;
    V.row(3) << 0.0, -0.5;

    Matrix2I R;
    R.row(0) << cos(theta), -sin(theta);
    R.row(1) << sin(theta), cos(theta);
    for (const auto& i : R.reshaped()) {
        CHECK(!empty(i));
    }

    Eigen::Matrix<filib::Interval, 4, 2> RV = V * R.transpose();
    RV.rowwise() += RowVector2I(x, y);
    for (const auto& i : R.reshaped()) {
        CHECK(!empty(i));
        CHECK(std::isfinite(i.INF));
        CHECK(std::isfinite(i.SUP));
    }
}

TEST_CASE("Interval norm", "[interval][matrix]")
{
    for (int i = 0; i < 100; i++) {
        Eigen::Vector3d v = Eigen::Vector3d::Random();
        VectorMax3I vi = v.cast<filib::Interval>();
        CHECK(in(v.norm(), norm(vi)));
    }

    Vector2I v(filib::Interval(-1, 1), filib::Interval(0, 0));
    CHECK(v.squaredNorm().INF == Catch::Approx(-1)); // This is wrong
    CHECK(v.squaredNorm().SUP == Catch::Approx(1));

    filib::Interval n = norm(v);
    CHECK(n.INF == 0);
    CHECK(n.SUP == Catch::Approx(1));
}

#endif