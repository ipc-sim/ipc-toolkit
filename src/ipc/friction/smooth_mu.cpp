#include "smooth_mu.hpp"

#include <ipc/friction/smooth_friction_mollifier.hpp>
#include <ipc/math/math.hpp>

#include <cassert>
#include <cmath>
#include <Eigen/Core>

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

std::pair<double, double> anisotropic_mu_eff_sqrt_mu0_t0_sq_plus_mu1_t1_sq(
    Eigen::ConstRef<Eigen::Vector2d> tau_dir,
    Eigen::ConstRef<Eigen::Vector2d> mu_s_aniso,
    Eigen::ConstRef<Eigen::Vector2d> mu_k_aniso)
{
    // Elliptical model (L2 projection):
    // mu_eff = sqrt((mu0*t0)^2 + (mu1*t1)^2) where t = tau_dir
    const double mu_s_eff = std::sqrt(
        Math<double>::sqr(mu_s_aniso[0] * tau_dir[0]) +
        Math<double>::sqr(mu_s_aniso[1] * tau_dir[1]));

    const double mu_k_eff = std::sqrt(
        Math<double>::sqr(mu_k_aniso[0] * tau_dir[0]) +
        Math<double>::sqr(mu_k_aniso[1] * tau_dir[1]));

    return std::make_pair(mu_s_eff, mu_k_eff);
}

Eigen::Vector2d anisotropic_mu_eff_dtau(
    Eigen::ConstRef<Eigen::Vector2d> tau,
    Eigen::ConstRef<Eigen::Vector2d> mu_aniso,
    const double mu_eff)
{
    constexpr double eps = 1e-10;
    const double tau_norm = tau.norm();

    // Edge cases: return zero vector if ||tau|| ~ 0 or mu_eff ~ 0
    if (tau_norm < eps || mu_eff < eps) {
        return Eigen::Vector2d::Zero();
    }

    // d(mu_eff)/d(tau_i) = tau_i * (mu_i^2 - mu_eff^2) / (mu_eff * ||tau||^2)
    const double tau_norm_sq = tau_norm * tau_norm;
    const double mu_eff_sq = mu_eff * mu_eff;
    Eigen::Vector2d result;
    result[0] = tau[0] * (Math<double>::sqr(mu_aniso[0]) - mu_eff_sq)
        / (mu_eff * tau_norm_sq);
    result[1] = tau[1] * (Math<double>::sqr(mu_aniso[1]) - mu_eff_sq)
        / (mu_eff * tau_norm_sq);

    // Ensure result is finite (handle numerical edge cases)
    if (!std::isfinite(result[0]) || !std::isfinite(result[1])) {
        return Eigen::Vector2d::Zero();
    }

    return result;
}

Eigen::Vector2d compute_tau_dir_from_tau_aniso(
    Eigen::ConstRef<Eigen::Vector2d> tau_aniso)
{
    constexpr double tiny = 1e-10;
    const double tau_aniso_norm = tau_aniso.norm();

    if (tau_aniso_norm < tiny) {
        return Eigen::Vector2d(1.0, 0.0); // Default direction
    } else {
        return tau_aniso / tau_aniso_norm;
    }
}

std::pair<double, double> compute_anisotropic_mu_eff_from_tau_aniso(
    Eigen::ConstRef<Eigen::Vector2d> tau_aniso,
    Eigen::ConstRef<Eigen::Vector2d> mu_s_aniso,
    Eigen::ConstRef<Eigen::Vector2d> mu_k_aniso,
    const double mu_s_isotropic,
    const double mu_k_isotropic,
    const bool no_mu)
{
    // Check if direction-dependent anisotropic friction is enabled
    const bool is_anisotropic = mu_s_aniso.squaredNorm() > 0;

    if (is_anisotropic) {
        // Direction-dependent friction: compute effective mu based on
        // tau_aniso direction (after mu_aniso scaling). This combines both
        // mechanisms: mu_aniso velocity scaling AND mu_s_aniso/mu_k_aniso
        // direction-dependent coefficients.
        const Eigen::Vector2d tau_dir =
            compute_tau_dir_from_tau_aniso(tau_aniso);

        const auto [mu_s_eff, mu_k_eff] =
            anisotropic_mu_eff_sqrt_mu0_t0_sq_plus_mu1_t1_sq(
                tau_dir, mu_s_aniso, mu_k_aniso);

        return std::make_pair(
            no_mu ? 1.0 : mu_s_eff, no_mu ? 1.0 : mu_k_eff);
    } else {
        // Isotropic friction: use scalar mu_s/mu_k (mu_aniso scaling already
        // applied to tau_aniso)
        return std::make_pair(
            no_mu ? 1.0 : mu_s_isotropic, no_mu ? 1.0 : mu_k_isotropic);
    }
}

std::pair<Eigen::Vector2d, Eigen::Vector2d>
compute_anisotropic_mu_eff_derivatives(
    Eigen::ConstRef<Eigen::Vector2d> tau_aniso,
    Eigen::ConstRef<Eigen::Vector2d> mu_s_aniso,
    Eigen::ConstRef<Eigen::Vector2d> mu_k_aniso,
    const double mu_s_eff,
    const double mu_k_eff)
{
    // Compute derivatives of effective mu w.r.t. tau_aniso
    const Eigen::Vector2d dmu_s_eff_dtau =
        anisotropic_mu_eff_dtau(tau_aniso, mu_s_aniso, mu_s_eff);
    const Eigen::Vector2d dmu_k_eff_dtau =
        anisotropic_mu_eff_dtau(tau_aniso, mu_k_aniso, mu_k_eff);

    return std::make_pair(dmu_s_eff_dtau, dmu_k_eff_dtau);
}

} // namespace ipc