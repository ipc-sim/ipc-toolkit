#pragma once

namespace ipc {

/// @brief Smooth friction mollifier function.
///
/// \f\[
///     f_0(y)= \begin{cases}
///         -\frac{y^3}{3\epsilon_v^2} + \frac{y^2}{\epsilon_v}
///             + \frac{\epsilon_v}{3}, & |y| < \epsilon_v
///             \newline
///         y, & |y| \geq \epsilon_v
///     \end{cases}
/// \f\]
///
/// @param y The tangential relative speed.
/// @param eps_v Velocity threshold below which static friction force is applied.
/// @return The value of the mollifier function at y.
double smooth_friction_f0(const double y, const double eps_v);

/// @brief The first derivative of the smooth friction mollifier.
///
/// \f\[
///     f_1(y) = f_0'(y) = \begin{cases}
///         -\frac{y^2}{\epsilon_v^2}+\frac{2 y}{\epsilon_v}, & |y| < \epsilon_v
///         \newline 1, & |y| \geq \epsilon_v
///     \end{cases}
/// \f\]
///
/// @param y The tangential relative speed.
/// @param eps_v Velocity threshold below which static friction force is applied.
/// @return The value of the derivative of the smooth friction mollifier at y.
double smooth_friction_f1(const double y, const double eps_v);

/// @brief The second derivative of the smooth friction mollifier.
///
/// \f\[
///     f_2(y) = f_0''(y) = \begin{cases}
///         -\frac{2 y}{\epsilon_v^2}+\frac{2}{\epsilon_v}, & |y| < \epsilon_v
///         \newline 0, & |y| \geq \epsilon_v
///     \end{cases}
/// \f\]
///
/// @param y The tangential relative speed.
/// @param eps_v Velocity threshold below which static friction force is applied.
/// @return The value of the second derivative of the smooth friction mollifier at y.
double smooth_friction_f2(const double y, const double eps_v);

/// @brief Compute the derivative of the smooth friction mollifier divided by y (\f$\frac{f_0'(y)}{y}\f$).
///
/// \f\[
///     \frac{f_1(y)}{y} = \begin{cases}
///         -\frac{y}{\epsilon_v^2}+\frac{2}{\epsilon_v}, & |y| < \epsilon_v
///         \newline \frac{1}{y}, & |y| \geq \epsilon_v
///     \end{cases}
/// \f\]
///
/// @param y The tangential relative speed.
/// @param eps_v Velocity threshold below which static friction force is applied.
/// @return The value of the derivative of smooth_friction_f0 divided by y.
double smooth_friction_f1_over_x(const double y, const double eps_v);

/// @brief The derivative of f1 times y minus f1 all divided by y cubed.
///
/// \f\[
///     \frac{f_1'(y) y - f_1(y)}{y^3} = \begin{cases}
///         -\frac{1}{y \epsilon_v^2}, & |y| < \epsilon_v \newline
///         -\frac{1}{y^3}, & |y| \geq \epsilon_v
///     \end{cases}
/// \f\]
///
/// @param y The tangential relative speed.
/// @param eps_v Velocity threshold below which static friction force is applied.
/// @return The derivative of f1 times y minus f1 all divided by y cubed.
double
smooth_friction_f2_x_minus_f1_over_x3(const double y, const double eps_v);

} // namespace ipc
