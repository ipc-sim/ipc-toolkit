#pragma once

namespace ipc {

// Fang and Li et al. [2023]:

// -- Normal Adhesion ----------------------------------------------------------
/// @defgroup normal_adhesion Normal Adhesion
/// @{

/// @brief The normal adhesion potential.
/// @param d distance
/// @param dhat_p distance of largest adhesion force (\f$\hat{d}_p\f$) where \f$0 < \hat{d}_p < \hat{d}_a\f$
/// @param dhat_a adhesion activation distance (\f$\hat{d}_a\f$)
/// @param a2 adjustable parameter relating to the maximum derivative of a (\f$a_2\f$)
/// @return The normal adhesion potential.
double normal_adhesion_potential(
    const double d, const double dhat_p, const double dhat_a, const double a2);

/// @brief The first derivative of the normal adhesion potential wrt d.
/// @param d distance
/// @param dhat_p distance of largest adhesion force (\f$\hat{d}_p\f$) where \f$0 < \hat{d}_p < \hat{d}_a\f$
/// @param dhat_a adhesion activation distance (\f$\hat{d}_a\f$)
/// @param a2 adjustable parameter relating to the maximum derivative of a (\f$a_2\f$)
/// @return The first derivative of the normal adhesion potential wrt d.
double normal_adhesion_potential_first_derivative(
    const double d, const double dhat_p, const double dhat_a, const double a2);

/// @brief The second derivative of the normal adhesion potential wrt d.
/// @param d distance
/// @param dhat_p distance of largest adhesion force (\f$\hat{d}_p\f$) where \f$0 < \hat{d}_p < \hat{d}_a\f$
/// @param dhat_a adhesion activation distance (\f$\hat{d}_a\f$)
/// @param a2 adjustable parameter relating to the maximum derivative of a (\f$a_2\f$)
/// @return The second derivative of the normal adhesion potential wrt d.
double normal_adhesion_potential_second_derivative(
    const double d, const double dhat_p, const double dhat_a, const double a2);

/// @brief The maximum normal adhesion force magnitude.
/// @param dhat_p distance of largest adhesion force (\f$\hat{d}_p\f$) where \f$0 < \hat{d}_p < \hat{d}_a\f$
/// @param dhat_a adhesion activation distance (\f$\hat{d}_a\f$)
/// @param a2 adjustable parameter relating to the maximum derivative of a (\f$a_2\f$)
/// @return The maximum normal adhesion force magnitude.
double max_normal_adhesion_force_magnitude(
    const double dhat_p, const double dhat_a, const double a2);

/// @}

// -- Tangential Adhesion ------------------------------------------------------
/// @defgroup tangential_adhesion Tangential Adhesion
/// @{

/// @brief The tangential adhesion mollifier function.
/// @param y The tangential relative speed.
/// @param eps_a Velocity threshold below which static adhesion force is applied.
/// @return The tangential adhesion mollifier function at y.
double tangential_adhesion_f0(const double y, const double eps_a);

/// @brief The first derivative of the tangential adhesion mollifier function.
/// @param y The tangential relative speed.
/// @param eps_a Velocity threshold below which static adhesion force is applied.
/// @return The first derivative of the tangential adhesion mollifier function at y.
double tangential_adhesion_f1(const double y, const double eps_a);

/// @brief The second derivative of the tangential adhesion mollifier function.
/// @param y The tangential relative speed.
/// @param eps_a Velocity threshold below which static adhesion force is applied.
/// @return The second derivative of the tangential adhesion mollifier function at y.
double tangential_adhesion_f2(const double y, const double eps_a);

/// @brief The first derivative of the tangential adhesion mollifier function divided by y.
/// @param y The tangential relative speed.
/// @param eps_a Velocity threshold below which static adhesion force is applied.
/// @return The first derivative of the tangential adhesion mollifier function divided by y.
double tangential_adhesion_f1_over_x(const double y, const double eps_a);

/// @brief The second derivative of the tangential adhesion mollifier function times y minus the first derivative all divided by y cubed.
/// @param y The tangential relative speed.
/// @param eps_a Velocity threshold below which static adhesion force is applied.
/// @return The second derivative of the tangential adhesion mollifier function times y minus the first derivative all divided by y cubed.
double
tangential_adhesion_f2_x_minus_f1_over_x3(const double y, const double eps_a);

// ~~ Smooth μ variants ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// NOTE: Here a0, a1, and a2 refer to the mollifier functions above.

/// @brief Compute the value of the ∫ μ(y) a₁(y) dy, where a₁ is the first derivative of the smooth tangential adhesion mollifier.
/// @note The `a0`/`a1` are unrelated to the `a0`/`a1` in the normal adhesion.
/// @param y The tangential relative speed.
/// @param mu_s Coefficient of static adhesion.
/// @param mu_k Coefficient of kinetic adhesion.
/// @param eps_a Velocity threshold below which static adhesion force is applied.
/// @return The value of the integral at y.
double smooth_mu_a0(
    const double y, const double mu_s, const double mu_k, const double eps_a);

/// @brief Compute the value of the μ(y) a₁(y), where a₁ is the first derivative of the smooth tangential adhesion mollifier.
/// @note The `a1` is unrelated to the `a1` in the normal adhesion.
/// @param y The tangential relative speed.
/// @param mu_s Coefficient of static adhesion.
/// @param mu_k Coefficient of kinetic adhesion.
/// @param eps_a Velocity threshold below which static adhesion force is applied.
/// @return The value of the product at y.
double smooth_mu_a1(
    const double y, const double mu_s, const double mu_k, const double eps_a);

/// @brief Compute the value of d/dy (μ(y) a₁(y)), where a₁ is the first derivative of the smooth tangential adhesion mollifier.
/// @note The `a1`/`a2` are unrelated to the `a1`/`a2` in the normal adhesion.
/// @param y The tangential relative speed.
/// @param mu_s Coefficient of static adhesion.
/// @param mu_k Coefficient of kinetic adhesion.
/// @param eps_a Velocity threshold below which static adhesion force is applied.
/// @return The value of the derivative at y.
double smooth_mu_a2(
    const double y, const double mu_s, const double mu_k, const double eps_a);

/// @brief Compute the value of the μ(y) a₁(y) / y, where a₁ is the first derivative of the smooth tangential adhesion mollifier.
/// @note The `x` in the function name refers to the parameter `y`.
/// @note The `a1` is unrelated to the `a1` in the normal adhesion.
/// @param y The tangential relative speed.
/// @param mu_s Coefficient of static adhesion.
/// @param mu_k Coefficient of kinetic adhesion.
/// @param eps_a Velocity threshold below which static adhesion force is applied.
/// @return The value of the product at y.
double smooth_mu_a1_over_x(
    const double y, const double mu_s, const double mu_k, const double eps_a);

/// @brief Compute the value of the [(d/dy μ(y) a₁(y)) ⋅ y - μ(y) a₁(y)] / y³, where a₁ and a₂ are the first and second derivatives of the smooth tangential adhesion mollifier.
/// @note The `x` in the function name refers to the parameter `y`.
/// @note The `a1`/`a2` are unrelated to the `a1`/`a2` in the normal adhesion.
/// @param y The tangential relative speed.
/// @param mu_s Coefficient of static adhesion.
/// @param mu_k Coefficient of kinetic adhesion.
/// @param eps_a Velocity threshold below which static adhesion force is applied.
/// @return The value of the expression at y.
double smooth_mu_a2_x_minus_mu_a1_over_x3(
    const double y, const double mu_s, const double mu_k, const double eps_a);

/// @}

} // namespace ipc
