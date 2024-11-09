#include "smooth_friction_mollifier.hpp"
#include <cassert>
#include <cmath>
#include <algorithm>
#include <optional>

namespace ipc {

// Helper function to select the friction coefficient based on blend type
double select_mu(
    const double s,
    const double epsv,
    const double mu1,
    const double mu2,
    const std::optional<BlendType>& blend_type)
{
    switch (blend_type.value_or(BlendType::TRANSITION)) {
    case BlendType::MIN:
        return std::min(mu1, mu2);
    case BlendType::MAX:
        return std::max(mu1, mu2);
    case BlendType::AVG:
        return 0.5 * (mu1 + mu2);
    case BlendType::TRANSITION:
        return (std::abs(s) < epsv) ? mu2 : mu1;
    default:
        assert(false);
        return 0.0;
    }
}

// Main mollifier function
double f0_SF(
    const double s,
    const double epsv,
    const std::optional<double>& static_mu,
    const std::optional<double>& kinetic_mu,
    const std::optional<BlendType>& blend_type)
{
    assert(epsv > 0);
    double mu = (static_mu && kinetic_mu)
        ? select_mu(s, epsv, *static_mu, *kinetic_mu, blend_type)
        : 1.0;

    if (std::abs(s) >= epsv) {
        return s * mu;
    }
    return (s * s * (-s / (3 * epsv) + 1) / epsv + epsv / 3) * mu;
}

// Derivative of mollifier function divided by s
double f1_SF_over_x(
    const double s,
    const double epsv,
    const std::optional<double>& static_mu,
    const std::optional<double>& kinetic_mu,
    const std::optional<BlendType>& blend_type)
{
    assert(epsv > 0);
    double mu = (static_mu && kinetic_mu)
        ? select_mu(s, epsv, *static_mu, *kinetic_mu, blend_type)
        : 1.0;

    if (std::abs(s) >= epsv) {
        return mu / s;
    }
    return (-s / epsv + 2) / epsv * mu;
}

// Derivative of f1 times s minus f1 all divided by s cubed
double df1_x_minus_f1_over_x3(
    const double s,
    const double epsv,
    const std::optional<double>& static_mu,
    const std::optional<double>& kinetic_mu,
    const std::optional<BlendType>& blend_type)
{
    assert(epsv > 0);
    double mu = (static_mu && kinetic_mu)
        ? select_mu(s, epsv, *static_mu, *kinetic_mu, blend_type)
        : 1.0;

    if (std::abs(s) >= epsv) {
        return -mu / (s * s * s);
    }
    return -mu / (s * epsv * epsv);
}

} // namespace ipc
