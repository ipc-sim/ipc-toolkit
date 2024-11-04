#pragma once

namespace ipc {

/// @brief Smooth friction mollifier function.
///
/// \f[
///     f_0(s)= \begin{cases}
///         -\frac{s^3}{3\epsilon_v^2} + \frac{s^2}{\epsilon_v}
///             + \frac{\epsilon_v}{3}, & |s| < \epsilon_v
///             \newline
///         s, & |s| \geq \epsilon_v
///     \end{cases}
/// \f]
///
/// @param s The tangential relative speed.
/// @param epsv Mollifier parameter \f$\epsilon_v\f$.
/// @return The value of the mollifier function at s.
double f0_SF(const double s, const double epsv);

/// @brief Compute the derivative of f0_SF divided by s (\f$\frac{f_0'(s)}{s}\f$).
///
/// \f[
///     f_1(s) = f_0'(s) = \begin{cases}
///         -\frac{s^2}{\epsilon_v^2}+\frac{2 s}{\epsilon_v}, & |s| < \epsilon_v
///         \newline 1, & |s| \geq \epsilon_v
///     \end{cases}
/// \f]
///
/// \f[
///     \frac{f_1(s)}{s} = \begin{cases}
///         -\frac{s}{\epsilon_v^2}+\frac{2}{\epsilon_v}, & |s| < \epsilon_v
///         \newline \frac{1}{s}, & |s| \geq \epsilon_v
///     \end{cases}
/// \f]
///
/// @param s The tangential relative speed.
/// @param epsv Mollifier parameter \f$\epsilon_v\f$.
/// @return The value of the derivative of f0_SF divided by s.
double f1_SF_over_x(const double s, const double epsv);

/// @brief The derivative of f1 times s minus f1 all divided by s cubed.
///
/// \f[
///     \frac{f_1'(s) s - f_1(s)}{s^3} = \begin{cases}
///         -\frac{1}{s \epsilon_v^2}, & |s| < \epsilon_v \newline
///         -\frac{1}{s^3}, & |s| \geq \epsilon_v
///     \end{cases}
/// \f]
///
/// @param s The tangential relative speed.
/// @param epsv Mollifier parameter \f$\epsilon_v\f$.
/// @return The derivative of f1 times s minus f1 all divided by s cubed.
double df1_x_minus_f1_over_x3(const double s, const double epsv);

// --------------------------------------------------------------------------
// New pairwise friction mollifier functions using blended mu coefficients.

/// @brief Blending function for pairwise friction coefficients.
/// @param mu1 Friction coefficient for material 1.
/// @param mu2 Friction coefficient for material 2.
/// @return The blended friction coefficient (default is the average of mu1 and mu2).
double blend_mu(const double mu1, const double mu2);

/// @brief Pairwise smooth friction mollifier function.
///
/// \f[
///     f_0(s)= \mu \begin{cases}
///         -\frac{s^3}{3\epsilon_v^2} + \frac{s^2}{\epsilon_v}
///             + \frac{\epsilon_v}{3}, & |s| < \epsilon_v
///             \newline
///         s, & |s| \geq \epsilon_v
///     \end{cases}
/// \f]
///
/// @param s The tangential relative speed.
/// @param epsv Mollifier parameter \f$\epsilon_v\f$.
/// @param mu1 Friction coefficient for material 1.
/// @param mu2 Friction coefficient for material 2.
/// @return The value of the mollifier function at s, using the blended mu.
double f0_SF_pairwise(
    const double s, const double epsv, const double mu1, const double mu2);

/// @brief Compute the derivative of f0_SF_pairwise divided by s (\f$\frac{f_0'(s)}{s}\f$).
///
/// \f[
///     f_1(s) = \mu f_0'(s) = \begin{cases}
///         -\frac{s^2}{\epsilon_v^2}+\frac{2 s}{\epsilon_v}, & |s| < \epsilon_v
///         \newline 1, & |s| \geq \epsilon_v
///     \end{cases}
/// \f]
///
/// \f[
///     \frac{f_1(s)}{s} = \mu \begin{cases}
///         -\frac{s}{\epsilon_v^2}+\frac{2}{\epsilon_v}, & |s| < \epsilon_v
///         \newline \frac{1}{s}, & |s| \geq \epsilon_v
///     \end{cases}
/// \f]
///
/// @param s The tangential relative speed.
/// @param epsv Mollifier parameter \f$\epsilon_v\f$.
/// @param mu1 Friction coefficient for material 1.
/// @param mu2 Friction coefficient for material 2.
/// @return The value of the derivative of f0_SF divided by s, using blended mu.
double f1_SF_over_x_pairwise(
    const double s, const double epsv, const double mu1, const double mu2);

/// @brief The derivative of f1 times s minus f1 all divided by s cubed, using blended mu.
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
/// @param mu1 Friction coefficient for material 1.
/// @param mu2 Friction coefficient for material 2.
/// @return The derivative of f1 times s minus f1 all divided by s cubed, using blended mu.
double df1_x_minus_f1_over_x3_pairwise(
    const double s, const double epsv, const double mu1, const double mu2);

// --------------------------------------------------------------------------
// Pairwise friction mollifier functions without blending, with kinetic/static transition.

/// @brief Pairwise smooth friction mollifier function without blending.
/// The function uses kinetic friction before transition (i.e., when \f$s < \epsilon_v\f$),
/// and static friction after the transition.
///
/// \f[
///     f_0(s)= \mu \begin{cases}
///         -\frac{s^3}{3\epsilon_v^2} + \frac{s^2}{\epsilon_v}
///             + \frac{\epsilon_v}{3}, & |s| < \epsilon_v
///             \newline
///         s, & |s| \geq \epsilon_v
///     \end{cases}
/// \f]
///
/// @param s The tangential relative speed.
/// @param epsv Mollifier parameter \f$\epsilon_v\f$.
/// @param static_mu Static friction coefficient.
/// @param kinetic_mu Kinetic friction coefficient.
/// @return The value of the mollifier function at s, using the correct friction coefficient.
double f0_SF_pairwise_transition(
    const double s, const double epsv, const double static_mu, const double kinetic_mu);

/// @brief Compute the derivative of f0_SF_pairwise_transition divided by s (\f$\frac{f_0'(s)}{s}\f$).
/// Uses kinetic friction before transition (i.e., when \f$s < \epsilon_v\f$),
/// and static friction after the transition.
///
/// @param s The tangential relative speed.
/// @param epsv Mollifier parameter \f$\epsilon_v\f$.
/// @param static_mu Static friction coefficient.
/// @param kinetic_mu Kinetic friction coefficient.
/// @return The value of the derivative of f0_SF divided by s, using the correct friction coefficient.
double f1_SF_over_x_pairwise_transition(
    const double s, const double epsv, const double static_mu, const double kinetic_mu);

/// @brief The derivative of f1 times s minus f1 all divided by s cubed,
/// using kinetic/static friction transition.
///
/// @param s The tangential relative speed.
/// @param epsv Mollifier parameter \f$\epsilon_v\f$.
/// @param static_mu Static friction coefficient.
/// @param kinetic_mu Kinetic friction coefficient.
/// @return The derivative of f1 times s minus f1 all divided by s cubed.
double df1_x_minus_f1_over_x3_pairwise_transition(
    const double s, const double epsv, const double static_mu, const double kinetic_mu);

} // namespace ipc
