#pragma once

namespace ipc {

/// @brief Smooth friction mollifier function.
///
/// \f\[
///     f_0(y)= \begin{cases}
///         -\frac{x^3}{3\epsilon_v^2 h^2} + \frac{x^2}{\epsilon_vh}
///             + \frac{\epsilon_v h}{3}, & |x| \in\left[0, \epsilon_v h\right)
///             \newline
///         x, & |x| \geq \epsilon_v h
///     \end{cases}
/// \f\]
///
/// @param x The tangential relative speed.
/// @param epsv_times_h Mollifier parameter \f$\epsilon_v h\f$.
/// @return The value of the mollifier function at x.
double f0_SF(const double x, const double epsv_times_h);

/// @brief Compute the derivative of f0_SF divided by x (\f$\frac{f_0'(x)}{x}\f$).
///
/// \f\[
///     f_1(x) = f_0'(x) = \begin{cases}
///         -\frac{x^2}{\epsilon_v^2 h^2}+\frac{2 x}{\epsilon_v h}, & |x|
///             \in\left[0, \epsilon_v h \right) \newline
///         1, & |x| \geq h \epsilon_v
///     \end{cases}
/// \f\]
///
/// \f\[
///     \frac{f_1(x)}{x} = \begin{cases}
///         -\frac{x}{\epsilon_v^2 h^2}+\frac{2}{\epsilon_v h}, & |x|
///             \in\left[0, \epsilon_v h \right) \newline
///         \frac{1}{x}, & |x| \geq h \epsilon_v
///     \end{cases}
/// \f\]
///
/// @param x The tangential relative speed.
/// @param epsv_times_h Mollifier parameter \f$\epsilon_v h\f$.
/// @return The value of the derivative of f0_SF divided by x.
double f1_SF_over_x(const double x, const double epsv_times_h);

/// @brief The derivative of f1 times x minus f1 all divided by x cubed.
///
/// \f\[
///     \frac{f_1'(x) x - f_1(x)}{x^3} = \begin{cases}
///         -\frac{1}{x \epsilon_v^2 h^2}, & |x|
///             \in\left[0, \epsilon_v h \right) \newline
///         -\frac{1}{x^3}, & |x| \geq h \epsilon_v
///     \end{cases}
/// \f\]
///
/// @param x The tangential relative speed.
/// @param epsv_times_h Mollifier parameter \f$\epsilon_v h\f$.
/// @return The derivative of f1 times x minus f1 all divided by x cubed.
double df1_x_minus_f1_over_x3(const double x, const double epsv_times_h);

} // namespace ipc
