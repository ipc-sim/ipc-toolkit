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
        const double z = std::abs(y) / eps_v;
        if (std::abs(y) < 0.5 * eps_v) {
            return 2 * (mu_k - mu_s) * z * z + mu_s;
        } else {
            return -2 * (mu_k - mu_s) * (z * (z - 2) + 1) + mu_k;
        }
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
        const double z = std::abs(y) / eps_v;
        if (std::abs(y) < 0.5 * eps_v) {
            return 4 * (mu_k - mu_s) * z / eps_v;
        } else {
            return -4 * (mu_k - mu_s) * (z - 1) / eps_v;
        }
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
        if (std::abs(y) < 0.5 * eps_v) {
            return y * z
                * (z * (z * (1 - 0.4 * z) * delta_mu - mu_s / 3.0) + mu_s)
                + (9.0 / 16.0) * eps_v * mu_k - (11.0 / 48.0) * eps_v * mu_s;
        } else {
            return y * z
                * (z
                       * (z * (0.4 * z - 2) * delta_mu
                          + (3 * mu_k - (10.0 / 3.0) * mu_s))
                   + (2 * mu_s - mu_k))
                + 0.6 * eps_v * mu_k - (4.0 / 15.0) * eps_v * mu_s;
        }
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

double smooth_mu_f2_x_minus_mu_f1_over_x3(
    const double y, const double mu_s, const double mu_k, const double eps_v)
{
    assert(eps_v > 0);
    if (mu_s == mu_k || std::abs(y) >= eps_v) {
        // If the static and kinetic friction coefficients are equal,
        // simplify.
        return mu_k * smooth_friction_f2_x_minus_f1_over_x3(y, eps_v);
    } else {
        const double delta_mu = mu_k - mu_s;
        const double z = 1 / eps_v;
        if (std::abs(y) < 0.5 * eps_v) {
            return z * z * (z * (8 - 6 * y * z) * delta_mu - mu_s / y);
        } else {
            return z * z
                * (z * (6 * y * z - 16) * delta_mu
                   + (9 * mu_k - 10 * mu_s) / y);
        }
    }
}

} // namespace ipc