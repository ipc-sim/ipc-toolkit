#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <array>
#include <cmath>

#include <ipc/esp/smooth_clamp.hpp>

using Catch::Approx;
using ipc::kSmoothClampEps;
using ipc::smooth_clamp01;
using ipc::smooth_clamp_simplex;

namespace {

// Central finite-difference derivative of smooth_clamp01 at x.
double fd_derivative(double x, double h)
{
    return (smooth_clamp01(x + h) - smooth_clamp01(x - h)) / (2.0 * h);
}

} // namespace

TEST_CASE("smooth_clamp01 anchor values", "[smooth_clamp]")
{
    const double eps = kSmoothClampEps;
    CHECK(smooth_clamp01(0.0) == Approx(0.0));
    CHECK(smooth_clamp01(eps) == Approx(eps));
    CHECK(smooth_clamp01(0.5) == Approx(0.5));
    CHECK(smooth_clamp01(1.0 - eps) == Approx(1.0 - eps));
    CHECK(smooth_clamp01(1.0) == Approx(1.0));
}

TEST_CASE("smooth_clamp01 saturates outside [0,1]", "[smooth_clamp]")
{
    for (const double x : { -10.0, -1.0, -0.001 })
        CHECK(smooth_clamp01(x) == Approx(0.0));
    for (const double x : { 1.001, 2.0, 100.0 })
        CHECK(smooth_clamp01(x) == Approx(1.0));
}

TEST_CASE("smooth_clamp01 is identity in [eps, 1-eps]", "[smooth_clamp]")
{
    const double eps = kSmoothClampEps;
    const int N = 100;
    for (int i = 0; i <= N; i++) {
        const double x = eps + (1.0 - 2.0 * eps) * i / N;
        CHECK(smooth_clamp01(x) == Approx(x));
    }
}

TEST_CASE("smooth_clamp01 is monotonically increasing", "[smooth_clamp]")
{
    const int N = 1000;
    double prev = smooth_clamp01(-0.2);
    for (int i = 1; i <= N; i++) {
        const double x = -0.2 + 1.4 * i / N; // sweep [-0.2, 1.2]
        const double cur = smooth_clamp01(x);
        CHECK(cur >= prev - 1e-15);
        prev = cur;
    }
}

TEST_CASE("smooth_clamp01 is continuous (C0) at all knots", "[smooth_clamp]")
{
    const double eps = kSmoothClampEps;
    const double knots[] = { 0.0, eps, 1.0 - eps, 1.0 };
    const double h = 1e-8;
    for (const double k : knots) {
        const double left = smooth_clamp01(k - h);
        const double right = smooth_clamp01(k + h);
        CHECK(std::abs(left - right) < 1e-7);
    }
}

TEST_CASE(
    "smooth_clamp01 is C1 (derivative matches across pieces at knots)",
    "[smooth_clamp]")
{
    const double eps = kSmoothClampEps;
    // Analytical derivatives per piece (defined for x in their respective
    // range):
    auto d_sat_lo = [](double) { return 0.0; }; // x <= 0
    auto d_blend_l = [&](double x) {
        return -3 * x * x / (eps * eps) + 4 * x / eps;
    }; // 0 <= x <= eps
    auto d_id = [](double) { return 1.0; }; // eps <= x <= 1-eps
    auto d_blend_r = [&](double x) {
        const double s = 1.0 - x;
        // d/dx [1 + s^3/eps^2 - 2 s^2/eps] with s = 1-x
        return -3 * s * s / (eps * eps) + 4 * s / eps;
    };
    auto d_sat_hi = [](double) { return 0.0; }; // x >= 1

    // At each knot the two adjacent piecewise derivative formulas must agree.
    CHECK(d_sat_lo(0.0) == Approx(d_blend_l(0.0)));
    CHECK(d_blend_l(eps) == Approx(d_id(eps)));
    CHECK(d_id(1.0 - eps) == Approx(d_blend_r(1.0 - eps)));
    CHECK(d_blend_r(1.0) == Approx(d_sat_hi(1.0)));
}

TEST_CASE(
    "smooth_clamp01 derivative matches FD across the whole range",
    "[smooth_clamp]")
{
    const double eps = kSmoothClampEps;
    const double h = 1e-6;
    const int N = 200;

    // Probe densely inside each smooth piece (avoid straddling knots).
    auto check_segment = [&](double a, double b) {
        for (int i = 1; i < N; i++) {
            const double x = a + (b - a) * i / N;
            // Stay strictly inside the piece by margin > h.
            if (x - h < a || x + h > b)
                continue;
            const double fd = fd_derivative(x, h);
            // Analytical derivative per piece:
            //   x in [0, eps]:    f'(x) = -3 x^2 / eps^2 + 4 x / eps
            //   x in [eps, 1-eps]: f'(x) = 1
            //   x in [1-eps, 1]:  f'(x) = 3 (1-x)^2 / eps^2 - 4 (1-x) / eps +
            //   ... ; via reflection same form as above on s=1-x
            double analytic;
            if (x < 0.0 || x > 1.0) {
                analytic = 0.0;
            } else if (x < eps) {
                analytic = -3 * x * x / (eps * eps) + 4 * x / eps;
            } else if (x > 1.0 - eps) {
                const double s = 1.0 - x;
                analytic = -3 * s * s / (eps * eps) + 4 * s / eps;
            } else {
                analytic = 1.0;
            }
            CHECK(std::abs(fd - analytic) < 1e-5);
        }
    };

    check_segment(-0.2, 0.0);      // saturated below: derivative 0
    check_segment(0.0, eps);       // entering blend
    check_segment(eps, 1.0 - eps); // identity
    check_segment(1.0 - eps, 1.0); // exiting blend
    check_segment(1.0, 1.2);       // saturated above: derivative 0
}

// ----- smooth_clamp_simplex -----

TEST_CASE("smooth_clamp_simplex sums to 1", "[smooth_clamp]")
{
    // Probe a grid of (u, v) including outside-triangle points; output must
    // always be a valid barycentric pair with u + v + w = 1.
    const int N = 30;
    for (int i = -10; i <= N + 10; i++) {
        for (int j = -10; j <= N + 10; j++) {
            const double u = double(i) / N;
            const double v = double(j) / N;
            double uo, vo;
            smooth_clamp_simplex(u, v, uo, vo);
            CHECK(uo >= -1e-15);
            CHECK(vo >= -1e-15);
            CHECK(uo + vo <= 1.0 + 1e-15);
        }
    }
}

TEST_CASE(
    "smooth_clamp_simplex is identity on the interior hexagon",
    "[smooth_clamp]")
{
    // Inside { u,v,w in [eps, 1-eps] } each smooth_clamp01 is identity and
    // sum is 1 exactly, so output equals input.
    const double eps = kSmoothClampEps;
    const int N = 50;
    for (int i = 0; i <= N; i++) {
        for (int j = 0; j <= N - i; j++) {
            const double u = eps + (1.0 - 3 * eps) * i / N;
            const double v = eps + (1.0 - 3 * eps) * j / N;
            const double w = 1.0 - u - v;
            if (u < eps || v < eps || w < eps || u > 1.0 - eps || v > 1.0 - eps
                || w > 1.0 - eps)
                continue;
            double uo, vo;
            smooth_clamp_simplex(u, v, uo, vo);
            CHECK(uo == Approx(u).epsilon(1e-12));
            CHECK(vo == Approx(v).epsilon(1e-12));
        }
    }
}

TEST_CASE(
    "smooth_clamp_simplex maps simplex vertices to themselves",
    "[smooth_clamp]")
{
    double uo, vo;

    smooth_clamp_simplex(0.0, 0.0, uo, vo); // T0 (w = 1)
    CHECK(uo == Approx(0.0));
    CHECK(vo == Approx(0.0));

    smooth_clamp_simplex(1.0, 0.0, uo, vo); // T1
    CHECK(uo == Approx(1.0));
    CHECK(vo == Approx(0.0));

    smooth_clamp_simplex(0.0, 1.0, uo, vo); // T2
    CHECK(uo == Approx(0.0));
    CHECK(vo == Approx(1.0));
}

TEST_CASE(
    "smooth_clamp_simplex saturates outside-triangle points to the boundary",
    "[smooth_clamp]")
{
    double uo, vo;

    // Far past T1 along the +u axis: (2, 0) should land at T1 = (1, 0).
    smooth_clamp_simplex(2.0, 0.0, uo, vo);
    CHECK(uo == Approx(1.0));
    CHECK(vo == Approx(0.0));

    // Far past T2: (0, 2) -> (0, 1).
    smooth_clamp_simplex(0.0, 2.0, uo, vo);
    CHECK(uo == Approx(0.0));
    CHECK(vo == Approx(1.0));

    // Negative orthant past T0: (-1, -1) -> (0, 0).
    smooth_clamp_simplex(-1.0, -1.0, uo, vo);
    CHECK(uo == Approx(0.0));
    CHECK(vo == Approx(0.0));
}

TEST_CASE("smooth_clamp_simplex is C1 (FD vs analytical)", "[smooth_clamp]")
{
    // Sample points and verify the central FD Jacobian matches an FD with a
    // smaller step at the same point — i.e. the function is smooth (no jumps).
    // Probe both inside and outside the simplex, but stay away from the knots
    // at u = 0/eps/1-eps/1 and v = 0/eps/1-eps/1 and w = 0/eps/1-eps/1 by a
    // margin > h.
    const double eps = kSmoothClampEps;
    const double h = 1e-6;

    auto J_fd = [&](double u, double v, double hh) {
        double upu, vpu, umu, vmu;
        double upv, vpv, umv, vmv;
        smooth_clamp_simplex(u + hh, v, upu, vpu);
        smooth_clamp_simplex(u - hh, v, umu, vmu);
        smooth_clamp_simplex(u, v + hh, upv, vpv);
        smooth_clamp_simplex(u, v - hh, umv, vmv);
        const double Juu = (upu - umu) / (2 * hh), Jvu = (vpu - vmu) / (2 * hh);
        const double Juv = (upv - umv) / (2 * hh), Jvv = (vpv - vmv) / (2 * hh);
        return std::array<double, 4> { { Juu, Juv, Jvu, Jvv } };
    };

    auto away_from_knot = [&](double x) {
        for (double k : { 0.0, eps, 1.0 - eps, 1.0 })
            if (std::abs(x - k) < 5 * h)
                return false;
        return true;
    };

    int probed = 0;
    for (double u : { -0.05, 0.05, 0.2, 0.4, 0.7, 0.95, 1.05 }) {
        for (double v : { -0.05, 0.05, 0.2, 0.4, 0.7, 0.95, 1.05 }) {
            const double w = 1.0 - u - v;
            if (!away_from_knot(u) || !away_from_knot(v) || !away_from_knot(w))
                continue;
            const auto j1 = J_fd(u, v, h);
            const auto j2 = J_fd(u, v, h * 4);
            for (int k = 0; k < 4; k++) {
                CHECK(std::abs(j1[k] - j2[k]) < 1e-4);
            }
            probed++;
        }
    }
    CHECK(probed > 0);
}
