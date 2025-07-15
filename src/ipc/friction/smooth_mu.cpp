#include "smooth_mu.hpp"

#include <ipc/friction/smooth_friction_mollifier.hpp>

#include <cassert>
#include <cmath>

namespace ipc {

double smooth_mu(
    const double y, const double mu_s, const double mu_k, const double eps_v)
{
    assert(eps_v > 0);
    if (mu_s == mu_k || std::abs(y) >= eps_v) {
        // If the static and kinetic friction coefficients are equal, simplify.
        return mu_k;
    } else {
        const double y_over_eps_v = std::abs(y) / eps_v;
        return (mu_k - mu_s) * (3 - 2 * y_over_eps_v) * y_over_eps_v
            * y_over_eps_v
            + mu_s;
    }
}

double smooth_mu_derivative(
    const double y, const double mu_s, const double mu_k, const double eps_v)
{
    assert(eps_v > 0);
    if (mu_s == mu_k || std::abs(y) >= eps_v) {
        // If the static and kinetic friction coefficients are equal, simplify.
        return 0;
    } else {
        const double y_over_eps_v = std::abs(y) / eps_v;
        return (mu_k - mu_s) * (6 - 6 * y_over_eps_v) * y_over_eps_v / eps_v;
    }
}

double smooth_mu_f0(
    const double y, const double mu_s, const double mu_k, const double eps_v)
{
    assert(eps_v > 0);
    if (mu_s == mu_k || std::abs(y) >= eps_v) {
        // If the static and kinetic friction coefficients are equal, simplify.
        return mu_k * smooth_friction_f0(y, eps_v);
    } else {
        const double delta_mu = mu_k - mu_s;
        const double z = std::abs(y) / eps_v;
        return y * z
            * (z * (z * (z * (z / 3 - 1.4) + 1.5) * delta_mu - mu_s / 3) + mu_s)
            + eps_v * (17 * mu_k - 7 * mu_s) / 30;
    }
}

double smooth_mu_f1(
    const double y, const double mu_s, const double mu_k, const double eps_v)
{
    // This is a known formulation: μ(y) f₁(y)
    return smooth_mu(y, mu_s, mu_k, eps_v) * smooth_friction_f1(y, eps_v);
}

double smooth_mu_f2(
    const double y, const double mu_s, const double mu_k, const double eps_v)
{
    // Apply the chain rule:
    return smooth_mu_derivative(y, mu_s, mu_k, eps_v)
        * smooth_friction_f1(y, eps_v)
        + smooth_mu(y, mu_s, mu_k, eps_v) * smooth_friction_f2(y, eps_v);
}

double smooth_mu_f1_over_x(
    const double y, const double mu_s, const double mu_k, const double eps_v)
{
    // This is a known formulation: μ(y) f₁(y) / y
    // where we use the robust division by y to avoid division by zero.
    return smooth_mu(y, mu_s, mu_k, eps_v)
        * smooth_friction_f1_over_x(y, eps_v);
}

double smooth_mu_f2_x_minus_f1_over_x3(
    const double y, const double mu_s, const double mu_k, const double eps_v)
{
    assert(eps_v > 0);
    if (mu_s == mu_k || std::abs(y) >= eps_v) {
        // If the static and kinetic friction coefficients are equal, simplify.
        return mu_k * smooth_friction_f2_x_minus_f1_over_x3(y, eps_v);
    } else {
        const double delta_mu = mu_k - mu_s;
        const double z = 1 / eps_v;
        return z * z
            * (z * (z * y * (z * y * 8 - 21) + 12) * delta_mu - mu_s / y);
    }
}

} // namespace ipc