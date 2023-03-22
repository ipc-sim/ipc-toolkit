#pragma once

namespace ipc {

/// @brief Smooth friction mollifier function.
///
/// \f\[
///     f_0(s)= \begin{cases}
///         -\frac{s^3}{3\epsilon_v^2} + \frac{s^2}{\epsilon_v}
///             + \frac{\epsilon_v}{3}, & |s| < \epsilon_v
///             \newline
///         s, & |s| \geq \epsilon_v
///     \end{cases}
/// \f\]
///
/// @param s The tangential relative speed.
/// @param epsv Mollifier parameter \f$\epsilon_v\f$.
/// @return The value of the mollifier function at s.
double f0_SF(const double s, const double epsv);

/// @brief Compute the derivative of f0_SF divided by s (\f$\frac{f_0'(s)}{s}\f$).
///
/// \f\[
///     f_1(s) = f_0'(s) = \begin{cases}
///         -\frac{s^2}{\epsilon_v^2}+\frac{2 s}{\epsilon_v}, & |s| < \epsilon_v
///         \newline 1, & |s| \geq \epsilon_v
///     \end{cases}
/// \f\]
///
/// \f\[
///     \frac{f_1(s)}{s} = \begin{cases}
///         -\frac{s}{\epsilon_v^2}+\frac{2}{\epsilon_v}, & |s| < \epsilon_v
///         \newline \frac{1}{s}, & |s| \geq \epsilon_v
///     \end{cases}
/// \f\]
///
/// @param s The tangential relative speed.
/// @param epsv Mollifier parameter \f$\epsilon_v\f$.
/// @return The value of the derivative of f0_SF divided by s.
double f1_SF_over_x(const double s, const double epsv);

/// @brief The derivative of f1 times s minus f1 all divided by s cubed.
///
/// \f\[
///     \frac{f_1'(s) s - f_1(s)}{s^3} = \begin{cases}
///         -\frac{1}{s \epsilon_v^2}, & |s| < \epsilon_v \newline
///         -\frac{1}{s^3}, & |s| \geq \epsilon_v
///     \end{cases}
/// \f\]
///
/// @param s The tangential relative speed.
/// @param epsv Mollifier parameter \f$\epsilon_v\f$.
/// @return The derivative of f1 times s minus f1 all divided by s cubed.
double df1_x_minus_f1_over_x3(const double s, const double epsv);

} // namespace ipc
