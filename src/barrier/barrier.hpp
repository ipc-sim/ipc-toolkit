/**
 * Barrier functions that grow to infinity as x -> 0+. Includes gradient and
 * hessian functions, too. These barrier functions can be used to impose
 * inequlity constraints on a function.
 */

#pragma once

#include <Eigen/Core>

namespace ipc {

/**
 * @brief Function that grows to infinity as x approaches 0 from the right.
 *
 * \f$b(d) = -(d-\hat{d})^2\ln(d / \hat{d})\f$
 *
 * @param x             The x value at which to evaluate.
 * @param s             Activation value of the barrier.
 * @param barrier_type  Barrier function type to use.
 * @return The value of the barrier function at x.
 */
template <typename T> T barrier(T d, double dhat);

double barrier_gradient(double d, double dhat);

double barrier_hessian(double d, double dhat);

} // namespace ipc

#include "barrier.tpp"
