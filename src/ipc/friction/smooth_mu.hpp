#pragma once

#include <ipc/utils/eigen_ext.hpp>

#include <utility> // for std::pair

namespace ipc {

/// @brief Smooth coefficient from static to kinetic friction.
/// @param y The tangential relative speed.
/// @param mu_s Coefficient of static friction.
/// @param mu_k Coefficient of kinetic friction.
/// @param eps_v Velocity threshold below which static friction force is applied.
/// @return The value of the μ at y.
double smooth_mu(
    const double y, const double mu_s, const double mu_k, const double eps_v);

/// @brief Compute the derivative of the smooth coefficient from static to kinetic friction.
/// @param y The tangential relative speed.
/// @param mu_s Coefficient of static friction.
/// @param mu_k Coefficient of kinetic friction.
/// @param eps_v Velocity threshold below which static friction force is applied.
/// @return The value of the derivative at y.
double smooth_mu_derivative(
    const double y, const double mu_s, const double mu_k, const double eps_v);

/// @brief Compute the value of the ∫ μ(y) f₁(y) dy, where f₁ is the first derivative of the smooth friction mollifier.
/// @param y The tangential relative speed.
/// @param mu_s Coefficient of static friction.
/// @param mu_k Coefficient of kinetic friction.
/// @param eps_v Velocity threshold below which static friction force is applied.
/// @return The value of the integral at y.
double smooth_mu_f0(
    const double y, const double mu_s, const double mu_k, const double eps_v);

/// @brief Compute the value of the μ(y) f₁(y), where f₁ is the first derivative of the smooth friction mollifier.
/// @param y The tangential relative speed.
/// @param mu_s Coefficient of static friction.
/// @param mu_k Coefficient of kinetic friction.
/// @param eps_v Velocity threshold below which static friction force is applied.
/// @return The value of the product at y.
double smooth_mu_f1(
    const double y, const double mu_s, const double mu_k, const double eps_v);

/// @brief Compute the value of d/dy (μ(y) f₁(y)), where f₁ is the first derivative of the smooth friction mollifier.
/// @param y The tangential relative speed.
/// @param mu_s Coefficient of static friction.
/// @param mu_k Coefficient of kinetic friction.
/// @param eps_v Velocity threshold below which static friction force is applied.
/// @return The value of the derivative at y.
double smooth_mu_f2(
    const double y, const double mu_s, const double mu_k, const double eps_v);

/// @brief Compute the value of the μ(y) f₁(y) / y, where f₁ is the first derivative of the smooth friction mollifier.
/// @note The `x` in the function name refers to the parameter `y`.
/// @param y The tangential relative speed.
/// @param mu_s Coefficient of static friction.
/// @param mu_k Coefficient of kinetic friction.
/// @param eps_v Velocity threshold below which static friction force is applied.
/// @return The value of the product at y.
double smooth_mu_f1_over_x(
    const double y, const double mu_s, const double mu_k, const double eps_v);

/// @brief Compute the value of the [(d/dy μ(y) f₁(y)) ⋅ y - μ(y) f₁(y)] / y³, where f₁ and f₂ are the first and second derivatives of the smooth friction mollifier.
/// @note The `x` in the function name refers to the parameter `y`.
/// @param y The tangential relative speed.
/// @param mu_s Coefficient of static friction.
/// @param mu_k Coefficient of kinetic friction.
/// @param eps_v Velocity threshold below which static friction force is applied.
/// @return The value of the expression at y.
double smooth_mu_f2_x_minus_mu_f1_over_x3(
    const double y, const double mu_s, const double mu_k, const double eps_v);

/// Elliptical L2 (matchstick cone) anisotropic friction. Call
/// anisotropic_x_from_tau_aniso, then anisotropic_mu_eff_f.
/// @see Erleben et al., CGF 2019, DOI 10.1111/cgf.13885;
///      https://github.com/erleben/matchstick

/// @brief Compute effective friction coefficients for elliptical anisotropy
///        (L2 projection): \f$\mu_{\text{eff}} = f(x) = \sqrt{(\mu_0 t_0)^2 +
///        (\mu_1 t_1)^2}\f$ at direction \f$x = \tau_{\text{dir}}\f$.
/// @note The \f$f\f$ in the name refers to the effective-μ formula; \f$x\f$ is the unit direction.
/// @details For anisotropic friction, the friction coefficient depends on the
///          direction of tangential velocity. This function computes the
///          effective friction coefficients along a given direction using the
///          elliptical L2 projection model:
///          \f$\mu_{\text{eff}} = \sqrt{(\mu_0 t_0)^2 + (\mu_1 t_1)^2}\f$,
///          where \f$t = \tau / \|\tau\|\f$ is the unit direction vector.
/// @param tau_dir Unit 2D direction of tangential velocity (tau / ||tau||).
///                Must be a unit vector.
/// @param mu_s_aniso Static friction ellipse axes (2D vector). Each component
///                   represents the friction coefficient along the
///                   corresponding tangent basis direction.
/// @param mu_k_aniso Kinetic friction ellipse axes (2D vector). Each component
///                   represents the friction coefficient along the
///                   corresponding tangent basis direction.
/// @return A pair containing (mu_s_eff, mu_k_eff), the effective static and
///         kinetic friction coefficients along the direction tau_dir.
/// @note If mu_s_aniso and mu_k_aniso are zero vectors, the function returns
///       (0, 0), which triggers compatible isotropic behavior.
/// @see anisotropic_x_from_tau_aniso, anisotropic_mu_eff_from_tau_aniso
[[nodiscard]] std::pair<double, double> anisotropic_mu_eff_f(
    Eigen::ConstRef<Eigen::Vector2d> tau_dir,
    Eigen::ConstRef<Eigen::Vector2d> mu_s_aniso,
    Eigen::ConstRef<Eigen::Vector2d> mu_k_aniso);

/// @brief Compute unit direction \f$x = \tau_{\text{aniso}} / \|\tau_{\text{aniso}}\|\f$ from tau_aniso, handling edge cases.
/// @param tau_aniso Anisotropically-scaled tangential velocity (2D vector).
/// @return Unit direction vector \f$x\f$. Returns (1, 0) if ||tau_aniso|| ≈ 0.
/// @see anisotropic_mu_eff_f, anisotropic_mu_eff_from_tau_aniso
[[nodiscard]] Eigen::Vector2d
anisotropic_x_from_tau_aniso(Eigen::ConstRef<Eigen::Vector2d> tau_aniso);

/// @brief Compute effective friction coefficients from tau_aniso, handling both
///        anisotropic and isotropic cases.
/// @details This function encapsulates the logic for determining whether to
///          use anisotropic or isotropic friction coefficients. If anisotropic
///          friction is enabled (mu_s_aniso.squaredNorm() > 0), it computes the
///          effective mu based on the direction of tau_aniso. Otherwise, it
///          returns the isotropic coefficients.
/// @param tau_aniso Anisotropically-scaled tangential velocity (2D vector).
/// @param mu_s_aniso Static friction ellipse axes (2D vector). Zero vector
///                   indicates isotropic friction.
/// @param mu_k_aniso Kinetic friction ellipse axes (2D vector). Zero vector
///                   indicates isotropic friction.
/// @param mu_s_isotropic Isotropic static friction coefficient (used when
///                       anisotropic is disabled).
/// @param mu_k_isotropic Isotropic kinetic friction coefficient (used when
///                       anisotropic is disabled).
/// @param no_mu If true, returns (1.0, 1.0) regardless of input coefficients.
/// @return A pair containing (mu_s, mu_k) to use in friction calculations.
/// @see anisotropic_x_from_tau_aniso, anisotropic_mu_eff_f
[[nodiscard]] std::pair<double, double> anisotropic_mu_eff_from_tau_aniso(
    Eigen::ConstRef<Eigen::Vector2d> tau_aniso,
    Eigen::ConstRef<Eigen::Vector2d> mu_s_aniso,
    Eigen::ConstRef<Eigen::Vector2d> mu_k_aniso,
    const double mu_s_isotropic,
    const double mu_k_isotropic,
    const bool no_mu = false);

} // namespace ipc