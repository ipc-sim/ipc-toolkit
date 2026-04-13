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

std::pair<double, double> anisotropic_mu_eff_f(
    Eigen::ConstRef<Eigen::Vector2d> tau_dir,
    Eigen::ConstRef<Eigen::Vector2d> mu_s_aniso,
    Eigen::ConstRef<Eigen::Vector2d> mu_k_aniso)
{
    // Elliptical model (L2 projection):
    // μ_eff = √((μ₀t₀)² + (μ₁t₁)²) where t = tau_dir
    const double mu_s_eff = mu_s_aniso.cwiseProduct(tau_dir).norm();
    const double mu_k_eff = mu_k_aniso.cwiseProduct(tau_dir).norm();
    return std::make_pair(mu_s_eff, mu_k_eff);
}

Eigen::Vector2d anisotropic_mu_eff_f_dtau(
    Eigen::ConstRef<Eigen::Vector2d> tau,
    Eigen::ConstRef<Eigen::Vector2d> mu_aniso,
    const double mu_eff)
{
    constexpr double EPS = 1e-10;
    const double tau_norm_sq = tau.squaredNorm();

    // Edge cases: return zero vector if ||tau|| ~ 0 or mu_eff ~ 0
    if (tau_norm_sq < EPS * EPS || mu_eff < EPS) {
        return Eigen::Vector2d::Zero();
    }

    // dμ_eff/dτᵢ = τᵢ * (μᵢ² - μ_eff²) / (μ_eff ‖τ‖²)
    const double mu_eff_sq = mu_eff * mu_eff;
    Eigen::Vector2d result = tau.array()
        * (mu_aniso.array().square() - mu_eff_sq) / (mu_eff * tau_norm_sq);

    // Ensure result is finite (handle numerical edge cases)
    return result.allFinite() ? result : Eigen::Vector2d::Zero();
}

Eigen::Vector2d
anisotropic_x_from_tau_aniso(Eigen::ConstRef<Eigen::Vector2d> tau_aniso)
{
    constexpr double EPS = 1e-10;
    const double tau_aniso_norm = tau_aniso.norm();

    if (tau_aniso_norm < EPS) {
        return Eigen::Vector2d(1.0, 0.0); // Default direction
    } else {
        return tau_aniso / tau_aniso_norm;
    }
}

std::pair<double, double> anisotropic_mu_eff_from_tau_aniso(
    Eigen::ConstRef<Eigen::Vector2d> tau_aniso,
    Eigen::ConstRef<Eigen::Vector2d> mu_s_aniso,
    Eigen::ConstRef<Eigen::Vector2d> mu_k_aniso,
    const double mu_s_isotropic,
    const double mu_k_isotropic,
    const bool no_mu)
{
    if (no_mu) {
        return std::make_pair(1.0, 1.0);
    }

    // Anisotropic when at least one of mu_s_aniso, mu_k_aniso is non-zero.
    // Compute mu_s_eff and mu_k_eff independently (ellipse or isotropic).
    const bool use_aniso_s = mu_s_aniso.squaredNorm() > 0;
    const bool use_aniso_k = mu_k_aniso.squaredNorm() > 0;

    if (!use_aniso_s && !use_aniso_k) {
        return std::make_pair(mu_s_isotropic, mu_k_isotropic);
    }

    const Eigen::Vector2d tau_dir = anisotropic_x_from_tau_aniso(tau_aniso);
    auto mu_eff = anisotropic_mu_eff_f(tau_dir, mu_s_aniso, mu_k_aniso);
    if (!use_aniso_s) {
        mu_eff.first = mu_s_isotropic;
    }
    if (!use_aniso_k) {
        mu_eff.second = mu_k_isotropic;
    }
    return mu_eff;
}

} // namespace ipc