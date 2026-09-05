#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <tests/simd_utils.hpp>

#ifdef IPC_TOOLKIT_WITH_SIMD

#include <ipc/adhesion/adhesion.hpp>

#include <array>
#include <string>

using namespace ipc;
using namespace ipc::tests;

TEST_CASE(
    "SIMD batch normal adhesion matches the scalar one lane-wise",
    "[adhesion][normal_adhesion][simd]")
{
    constexpr double DHAT_P = 1e-3;
    constexpr double DHAT_A = 2e-3;
    const double max_slope = GENERATE(-1.0, -1e3);

    // Distances landing in each piece: the quadratic below d̂ₚ, the second
    // quadratic between d̂ₚ and d̂ₐ, and the inactive region past d̂ₐ -- plus the
    // two breakpoints themselves, where the pieces must agree.
    constexpr std::array<double, 7> DS = { 0.0,    0.5e-3, DHAT_P, 1.5e-3,
                                           DHAT_A, 3e-3,   1.0 };

    auto check = [&](const std::string& name, auto&& f) {
        check_swept_lanes(
            name, DS,
            [&](const double d) { return f(d, DHAT_P, DHAT_A, max_slope); },
            [&](const Batch& d) {
                return f(d, Batch(DHAT_P), Batch(DHAT_A), Batch(max_slope));
            });
    };

    check("potential", [](auto d, auto dhat_p, auto dhat_a, auto a2) {
        return normal_adhesion_potential(d, dhat_p, dhat_a, a2);
    });
    check("first_derivative", [](auto d, auto dhat_p, auto dhat_a, auto a2) {
        return normal_adhesion_potential_first_derivative(
            d, dhat_p, dhat_a, a2);
    });
    check("second_derivative", [](auto d, auto dhat_p, auto dhat_a, auto a2) {
        return normal_adhesion_potential_second_derivative(
            d, dhat_p, dhat_a, a2);
    });
}

TEST_CASE(
    "SIMD batch tangential adhesion matches the scalar one lane-wise",
    "[adhesion][tangential_adhesion][simd]")
{
    const double eps_a = GENERATE(1e-3, 0.1, 1.0);

    // Speeds as multiples of eps_a. These functions clamp at y <= 0 rather than
    // mirroring on |y|, so the negative entry checks that clamp; zero is where
    // the `1/y` branch is singular, which only a batch evaluates.
    constexpr std::array<double, 8> MULTIPLES = { -1.0, 0.0,  0.25, 0.49,
                                                  0.5,  0.75, 1.0,  2.5 };
    std::array<double, MULTIPLES.size()> ys {};
    for (size_t i = 0; i < ys.size(); ++i) {
        ys[i] = MULTIPLES[i] * eps_a;
    }

    auto check = [&](const std::string& name, auto&& f) {
        check_swept_lanes(
            name, ys, [&](const double y) { return f(y, eps_a); },
            [&](const Batch& y) { return f(y, Batch(eps_a)); });
    };

    check("f0", [](auto y, auto e) { return tangential_adhesion_f0(y, e); });
    check("f1", [](auto y, auto e) { return tangential_adhesion_f1(y, e); });
    check("f2", [](auto y, auto e) { return tangential_adhesion_f2(y, e); });
    check("f1_over_x", [](auto y, auto e) {
        return tangential_adhesion_f1_over_x(y, e);
    });
}

TEST_CASE(
    "SIMD batch smooth mu adhesion variants match the scalar ones lane-wise",
    "[adhesion][smooth_mu][simd]")
{
    const double eps_a = GENERATE(1e-3, 0.1, 1.0);

    // The equal-coefficient pair is the fast path these functions short-circuit
    // on, and it has to agree with the general formulas it skips.
    const auto mus = GENERATE(
        std::pair { 0.5, 0.5 }, std::pair { 0.5, 0.1 }, std::pair { 0.1, 0.5 });
    const double mu_s = mus.first, mu_k = mus.second;

    // Non-negative only: smooth_mu_a2_x_minus_mu_a1_over_x3 asserts y >= 0.
    constexpr std::array<double, 7> MULTIPLES = { 0.0,  0.25, 0.49, 0.5,
                                                  0.75, 1.0,  2.5 };
    std::array<double, MULTIPLES.size()> ys {};
    for (size_t i = 0; i < ys.size(); ++i) {
        ys[i] = MULTIPLES[i] * eps_a;
    }

    auto check = [&](const std::string& name, auto&& f) {
        check_swept_lanes(
            name, ys, [&](const double y) { return f(y, mu_s, mu_k, eps_a); },
            [&](const Batch& y) {
                return f(y, Batch(mu_s), Batch(mu_k), Batch(eps_a));
            });
    };

    check("a0", [](auto y, auto s, auto k, auto e) {
        return smooth_mu_a0(y, s, k, e);
    });
    check("a1", [](auto y, auto s, auto k, auto e) {
        return smooth_mu_a1(y, s, k, e);
    });
    check("a2", [](auto y, auto s, auto k, auto e) {
        return smooth_mu_a2(y, s, k, e);
    });
    check("a1_over_x", [](auto y, auto s, auto k, auto e) {
        return smooth_mu_a1_over_x(y, s, k, e);
    });
    check("a2_x_minus_mu_a1_over_x3", [](auto y, auto s, auto k, auto e) {
        return smooth_mu_a2_x_minus_mu_a1_over_x3(y, s, k, e);
    });
}

#endif
