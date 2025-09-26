#pragma once

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

} // namespace ipc