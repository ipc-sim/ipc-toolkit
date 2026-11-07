#include "smooth_mu.hpp"

namespace ipc {

// The anisotropic helpers stay `double`-only for now, unlike the scalar μ
// functions in the header, because no batch caller exists for them yet. Two of
// them do branch per collision: `anisotropic_x_from_tau_aniso` guards a
// vanishing tangential velocity and `anisotropic_mu_eff_f` is a function of
// its direction. A batch friction path with anisotropic μ would therefore need
// to template those two as well. Only `anisotropic_mu_eff_from_tau_aniso`
// tests the material alone -- whether an ellipse axis is set, whether the
// caller disabled μ -- and can stay scalar regardless.

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
