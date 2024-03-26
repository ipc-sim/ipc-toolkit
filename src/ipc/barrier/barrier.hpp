// Barrier functions that grow to infinity as x -> 0+. Includes gradient and
// hessian functions, too. These barrier functions can be used to impose
// inequality constraints on a function.

#pragma once

namespace ipc {

/// Base class for barrier functions.
class Barrier {
public:
    Barrier() = default;
    virtual ~Barrier() = default;

    virtual double operator()(const double d, const double dhat) const = 0;
    virtual double
    first_derivative(const double d, const double dhat) const = 0;
    virtual double
    second_derivative(const double d, const double dhat) const = 0;
};

// ============================================================================
// Barrier functions from [Li et al. 2020]
// ============================================================================

/// @brief Function that grows to infinity as d approaches 0 from the right.
///
/// \f\[
///     b(d) = -(d-\hat{d})^2\ln\left(\frac{d}{\hat{d}}\right)
/// \f\]
///
/// @param d The distance.
/// @param dhat Activation distance of the barrier.
/// @return The value of the barrier function at d.
double barrier(const double d, const double dhat);

/// @brief Derivative of the barrier function.
///
/// \f\[
///     b'(d) = (\hat{d}-d) \left( 2\ln\left( \frac{d}{\hat{d}} \right) -
///     \frac{\hat{d}}{d} + 1\right)
/// \f\]
///
/// @param d The distance.
/// @param dhat Activation distance of the barrier.
/// @return The derivative of the barrier wrt d.
double barrier_first_derivative(const double d, const double dhat);

/// @brief Second derivative of the barrier function.
///
/// \f\[
///     b''(d) = \left( \frac{\hat{d}}{d} + 2 \right) \frac{\hat{d}}{d} -
///     2\ln\left( \frac{d}{\hat{d}} \right) - 3
/// \f\]
///
/// @param d The distance.
/// @param dhat Activation distance of the barrier.
/// @return The second derivative of the barrier wrt d.
double barrier_second_derivative(const double d, const double dhat);

/// @brief Smoothly clamped log barrier functions from [Li et al. 2020].
class ClampedLogBarrier : public Barrier {
public:
    ClampedLogBarrier() = default;

    /// @brief Function that grows to infinity as d approaches 0 from the right.
    ///
    /// \f\[
    ///     b(d) = -(d-\hat{d})^2\ln\left(\frac{d}{\hat{d}}\right)
    /// \f\]
    ///
    /// @param d The distance.
    /// @param dhat Activation distance of the barrier.
    /// @return The value of the barrier function at d.
    double operator()(const double d, const double dhat) const override
    {
        return barrier(d, dhat);
    }

    /// @brief Derivative of the barrier function.
    ///
    /// \f\[
    ///     b'(d) = (\hat{d}-d) \left( 2\ln\left( \frac{d}{\hat{d}} \right) -
    ///     \frac{\hat{d}}{d} + 1\right)
    /// \f\]
    ///
    /// @param d The distance.
    /// @param dhat Activation distance of the barrier.
    /// @return The derivative of the barrier wrt d.
    double first_derivative(const double d, const double dhat) const override
    {
        return barrier_first_derivative(d, dhat);
    }

    /// @brief Second derivative of the barrier function.
    ///
    /// \f\[
    ///     b''(d) = \left( \frac{\hat{d}}{d} + 2 \right) \frac{\hat{d}}{d} -
    ///     2\ln\left( \frac{d}{\hat{d}} \right) - 3
    /// \f\]
    ///
    /// @param d The distance.
    /// @param dhat Activation distance of the barrier.
    /// @return The second derivative of the barrier wrt d.
    double second_derivative(const double d, const double dhat) const override
    {
        return barrier_second_derivative(d, dhat);
    }
};

} // namespace ipc
