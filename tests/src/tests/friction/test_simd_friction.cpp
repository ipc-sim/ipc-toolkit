#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <tests/simd_utils.hpp>

#ifdef IPC_TOOLKIT_WITH_SIMD

#include <ipc/friction/smooth_friction_mollifier.hpp>
#include <ipc/friction/smooth_mu.hpp>

#include <array>
#include <string>

using namespace ipc;
using namespace ipc::tests;

namespace {

/// @brief Speeds placing a lane in each piece of the mollifier, as multiples of
/// eps_v: well inside, either side of the half-threshold the smooth-μ formulas
/// split on, just inside, exactly at, and past the threshold -- plus the
/// negative mirror, since these functions branch on |y| but the value itself
/// carries the sign. Zero is in the list because that is where the `1/y`
/// branches are singular: the scalar path never takes them there, while a batch
/// evaluates them anyway and must blend the infinity away.
constexpr std::array<double, 11> Y_MULTIPLES = { -2.0, -1.0, -0.75, -0.25,
                                                 0.0,  0.25, 0.49,  0.5,
                                                 0.75, 1.0,  2.0 };

} // namespace

TEST_CASE(
    "SIMD batch smooth friction mollifier matches the scalar one lane-wise",
    "[friction][mollifier][simd]")
{
    const double eps_v = GENERATE(1e-3, 0.1, 1.0);
    const auto ys = scaled(Y_MULTIPLES, eps_v);

    auto check = [&](const std::string& name, auto&& f) {
        check_swept_lanes_with(name, ys, f, eps_v);
    };

    check("f0", [](auto y, auto e) { return smooth_friction_f0(y, e); });
    check("f1", [](auto y, auto e) { return smooth_friction_f1(y, e); });
    check("f2", [](auto y, auto e) { return smooth_friction_f2(y, e); });
    check("f1_over_x", [](auto y, auto e) {
        return smooth_friction_f1_over_x(y, e);
    });
    check("f2_x_minus_f1_over_x3", [](auto y, auto e) {
        return smooth_friction_f2_x_minus_f1_over_x3(y, e);
    });
}

TEST_CASE(
    "SIMD batch smooth mu matches the scalar one lane-wise",
    "[friction][smooth_mu][simd]")
{
    const double eps_v = GENERATE(1e-3, 0.1, 1.0);

    // The equal-coefficient pair is the fast path these functions short-circuit
    // on, and it has to agree with the general formulas it skips.
    const auto mus = GENERATE(
        std::pair { 0.5, 0.5 }, std::pair { 0.5, 0.1 }, std::pair { 0.1, 0.5 });
    const double mu_s = mus.first, mu_k = mus.second;

    const auto ys = scaled(Y_MULTIPLES, eps_v);

    auto check = [&](const std::string& name, auto&& f) {
        check_swept_lanes_with(name, ys, f, mu_s, mu_k, eps_v);
    };

    check("mu", [](auto y, auto s, auto k, auto e) {
        return smooth_mu(y, s, k, e);
    });
    check("mu_derivative", [](auto y, auto s, auto k, auto e) {
        return smooth_mu_derivative(y, s, k, e);
    });
    check("mu_f0", [](auto y, auto s, auto k, auto e) {
        return smooth_mu_f0(y, s, k, e);
    });
    check("mu_f1", [](auto y, auto s, auto k, auto e) {
        return smooth_mu_f1(y, s, k, e);
    });
    check("mu_f2", [](auto y, auto s, auto k, auto e) {
        return smooth_mu_f2(y, s, k, e);
    });
    check("mu_f1_over_x", [](auto y, auto s, auto k, auto e) {
        return smooth_mu_f1_over_x(y, s, k, e);
    });
    check("mu_f2_x_minus_mu_f1_over_x3", [](auto y, auto s, auto k, auto e) {
        return smooth_mu_f2_x_minus_mu_f1_over_x3(y, s, k, e);
    });
}

#endif
