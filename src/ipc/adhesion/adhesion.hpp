#pragma once

namespace ipc {

// Fang and Li et al. [2023]:

// -- Normal Adhesion ----------------------------------------------------------

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

// -- Tangential Adhesion ------------------------------------------------------

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

} // namespace ipc
