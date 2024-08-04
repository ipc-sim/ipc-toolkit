#pragma once

namespace ipc {

// Fang and Li et al. [2023]:

// -- Normal Adhesion ----------------------------------------------------------

/// @brief The normal adhesion potential.
/// @param d distance
/// @param dhat_p distance of largest adhesion force (\f(\hat{d}_p\f)) (\f(0 < \hat{d}_p < \hat{d}_a\f))
/// @param dhat_a adhesion activation distance (\f(\hat{d}_a\f))
/// @param a2 adjustable parameter relating to the maximum derivative of a (\f(a_2\f))
/// @return The normal adhesion potential.
double normal_adhesion_potential(
    const double d, const double dhat_p, const double dhat_a, const double a2);

/// @brief The first derivative of the normal adhesion potential wrt d.
/// @param d distance
/// @param dhat_p distance of largest adhesion force (\f(\hat{d}_p\f)) (\f(0 < \hat{d}_p < \hat{d}_a\f))
/// @param dhat_a adhesion activation distance (\f(\hat{d}_a\f))
/// @param a2 adjustable parameter relating to the maximum derivative of a (\f(a_2\f))
/// @return The first derivative of the normal adhesion potential wrt d.
double normal_adhesion_potential_first_derivative(
    const double d, const double dhat_p, const double dhat_a, const double a2);

/// @brief The second derivative of the normal adhesion potential wrt d.
/// @param d distance
/// @param dhat_p distance of largest adhesion force (\f(\hat{d}_p\f)) (\f(0 < \hat{d}_p < \hat{d}_a\f))
/// @param dhat_a adhesion activation distance (\f(\hat{d}_a\f))
/// @param a2 adjustable parameter relating to the maximum derivative of a (\f(a_2\f))
/// @return The second derivative of the normal adhesion potential wrt d.
double normal_adhesion_potential_second_derivative(
    const double d, const double dhat_p, const double dhat_a, const double a2);

// -- Tangential Adhesion ------------------------------------------------------

/// @brief The tangent adhesion potential.
/// @param y
/// @param eps_a
/// @return
double f0_t(const double y, const double eps_a);

/// @brief The tangent adhesion potential gradient.
double f1_t(const double y, const double eps_a);

double df1_t(const double y, const double eps_a);

double f1_t_over_x(const double y, const double eps_a);

double df1_t_x_minus_f1_t_over_x3(const double y, const double eps_a);

} // namespace ipc
