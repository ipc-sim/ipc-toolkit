// Barrier functions that grow to infinity as x -> 0+. Includes gradient and
// hessian functions, too. These barrier functions can be used to impose
// inequlity constraints on a function.

#pragma once

#include <Eigen/Core>

namespace ipc {

/// @brief Function that grows to infinity as x approaches 0 from the right.
///
/// \f$b(d) = -(d-\hat{d})^2\ln(d / \hat{d})\f$
///
/// @param d The distance.
/// @param dhat Activation distance of the barrier.
/// @return The value of the barrier function at d.
template <typename T> T barrier(const T& d, double dhat);

/// @brief Derivative of the barrier function.
/// @param d The distance.
/// @param dhat Activation distance of the barrier.
/// @return The derivative of the barrier wrt d.
double barrier_gradient(double d, double dhat);

/// @brief Second derivative of the barrier function.
/// @param d The distance.
/// @param dhat Activation distance of the barrier.
/// @return The second derivative of the barrier wrt d.
double barrier_hessian(double d, double dhat);

} // namespace ipc

#include "barrier.tpp"
