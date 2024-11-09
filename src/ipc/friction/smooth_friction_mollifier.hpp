#pragma once
#include <optional>

namespace ipc {

enum class BlendType {
    AVG = 0,       // Average
    MIN = 1,       // Minimum
    MAX = 2,       // Maximum
    TRANSITION = 3 // Kinetic/static friction transition
};

/// @brief Smooth friction mollifier function for friction forces.
/// Applies friction blending or transition-based coefficients.
///
/// \f[
///     f_0(s) = \mu \begin{cases}
///         -\frac{s^3}{3\epsilon_v^2} + \frac{s^2}{\epsilon_v} + \frac{\epsilon_v}{3}, & |s| < \epsilon_v \newline
///         s, & |s| \geq \epsilon_v
///     \end{cases}
/// \f]
///
/// @param s The tangential relative speed.
/// @param epsv Mollifier parameter \f$\epsilon_v\f$.
/// @param static_mu Optional static friction coefficient.
/// @param kinetic_mu Optional kinetic friction coefficient.
/// @param blend_type Optional type of blending for coefficients.
/// @return The mollified friction force at relative speed s.
double f0_SF(
    const double s,
    const double epsv,
    const std::optional<double>& static_mu = std::nullopt,
    const std::optional<double>& kinetic_mu = std::nullopt,
    const std::optional<BlendType>& blend_type = std::nullopt);

/// @brief Derivative of f0_SF divided by s.
///
/// \f[
///     f_1(s) = f_0'(s) = \mu \begin{cases}
///         -\frac{s^2}{\epsilon_v^2} + \frac{2 s}{\epsilon_v}, & |s| < \epsilon_v \newline
///         1, & |s| \geq \epsilon_v
///     \end{cases}
/// \f]
///
/// \f[
///     \frac{f_1(s)}{s} = \mu \begin{cases}
///         -\frac{s}{\epsilon_v^2} + \frac{2}{\epsilon_v}, & |s| < \epsilon_v \newline
///         \frac{1}{s}, & |s| \geq \epsilon_v
///     \end{cases}
/// \f]
///
/// @param s The tangential relative speed.
/// @param epsv Mollifier parameter \f$\epsilon_v\f$.
/// @param static_mu Optional static friction coefficient.
/// @param kinetic_mu Optional kinetic friction coefficient.
/// @param blend_type Optional type of blending for coefficients.
/// @return The derivative of f0_SF with respect to s.
double f1_SF_over_x(
    const double s,
    const double epsv,
    const std::optional<double>& static_mu = std::nullopt,
    const std::optional<double>& kinetic_mu = std::nullopt,
    const std::optional<BlendType>& blend_type = std::nullopt);

/// @brief The derivative of f1 times s minus f1 all divided by s cubed.
///
/// \f[
///     \frac{f_1'(s) s - f_1(s)}{s^3} = \mu \begin{cases}
///         -\frac{1}{s \epsilon_v^2}, & |s| < \epsilon_v \newline
///         -\frac{1}{s^3}, & |s| \geq \epsilon_v
///     \end{cases}
/// \f]
///
/// @param s The tangential relative speed.
/// @param epsv Mollifier parameter \f$\epsilon_v\f$.
/// @param static_mu Optional static friction coefficient.
/// @param kinetic_mu Optional kinetic friction coefficient.
/// @param blend_type Optional type of blending for coefficients.
/// @return The derivative of f1 times s minus f1 all divided by s cubed.
double df1_x_minus_f1_over_x3(
    const double s,
    const double epsv,
    const std::optional<double>& static_mu = std::nullopt,
    const std::optional<double>& kinetic_mu = std::nullopt,
    const std::optional<BlendType>& blend_type = std::nullopt);

/// @brief Blending function for pairwise friction coefficients.
///
/// Uses blending to return a combined coefficient based on the type of blend:
///
/// \f[
///     \mu = \begin{cases}
///         \frac{\mu_1 + \mu_2}{2}, & \text{if type = AVG} \newline
///         \min(\mu_1, \mu_2), & \text{if type = MIN} \newline
///         \max(\mu_1, \mu_2), & \text{if type = MAX} \newline
///         \text{(depends on } s \text{ for kinetic/static transition)}, & \text{if type = TRANSITION}
///     \end{cases}
/// \f]
///
/// @param mu1 Friction coefficient for material 1.
/// @param mu2 Friction coefficient for material 2.
/// @param type The type of blending to use.
/// @return The blended friction coefficient (default is the average of mu1 and mu2).
double blend_mu(
    const double mu1,
    const double mu2,
    const std::optional<BlendType> type = std::nullopt);

} // namespace ipc
