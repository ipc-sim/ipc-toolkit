// Barrier functions that grow to infinity as x -> 0+. Includes gradient and
// hessian functions, too. These barrier functions can be used to impose
// inequality constraints on a function.

#pragma once

#include <ipc/config.hpp>

namespace ipc {

/// Base class for barrier functions.
template <typename T = double> class BarrierBase {
protected:
    using value_type = T;

public:
    BarrierBase() = default;
    virtual ~BarrierBase() = default;

    /// @brief Evaluate the barrier function.
    /// @param d Distance.
    /// @param dhat Activation distance of the barrier.
    /// @return The value of the barrier function at d.
    virtual T operator()(const T d, const T dhat) const = 0;

    /// @brief Evaluate the first derivative of the barrier function wrt d.
    /// @param d Distance.
    /// @param dhat Activation distance of the barrier.
    /// @return The value of the first derivative of the barrier function at d.
    virtual T first_derivative(const T d, const T dhat) const = 0;

    /// @brief Evaluate the second derivative of the barrier function wrt d.
    /// @param d Distance.
    /// @param dhat Activation distance of the barrier.
    /// @return The value of the second derivative of the barrier function at d.
    virtual T second_derivative(const T d, const T dhat) const = 0;

    /// @brief Get the units of the barrier function.
    /// Essentially, barrier(d, d̂) / units(d̂) should be dimensionless.
    /// @param dhat The activation distance of the barrier.
    /// @return The units of the barrier function.
    virtual T units(const T dhat) const = 0;
};

/// @brief Default barrier type.
using Barrier = BarrierBase<>;

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
template <typename T = double> T barrier(const T d, const T dhat);

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
template <typename T = double>
T barrier_first_derivative(const T d, const T dhat);

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
template <typename T = double>
T barrier_second_derivative(const T d, const T dhat);

/// @brief Smoothly clamped log barrier functions from [Li et al. 2020].
template <typename T = double> class ClampedLogBarrier : public BarrierBase<T> {
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
    T operator()(const T d, const T dhat) const override
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
    T first_derivative(const T d, const T dhat) const override
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
    T second_derivative(const T d, const T dhat) const override
    {
        return barrier_second_derivative(d, dhat);
    }

    /// @brief Get the units of the barrier function.
    /// @param dhat The activation distance of the barrier.
    /// @return The units of the barrier function.
    T units(const T dhat) const override
    {
        // (d - d̂)² = d̂² (d/d̂ - 1)²
        return dhat * dhat;
    }
};

// ============================================================================
// Normalized Barrier functions from [Li et al. 2023]
// ============================================================================

/// @brief Normalized barrier function from [Li et al. 2023].
template <typename BarrierT> class NormalizedBarrier : public BarrierT {
    using T = typename BarrierT::value_type;

public:
    NormalizedBarrier() = default;

    /// @brief Function that grows to infinity as d approaches 0 from the right.
    ///
    /// \f\[
    ///     b(d) =
    ///     -\left(\frac{d}{\hat{d}}-1\right)^2\ln\left(\frac{d}{\hat{d}}\right)
    /// \f\]
    ///
    /// @param d The distance.
    /// @param dhat Activation distance of the barrier.
    /// @return The value of the barrier function at d.
    T operator()(const T d, const T dhat) const override
    {
        return BarrierT::operator()(d / dhat, T(1));
    }

    /// @brief Derivative of the barrier function.
    ///
    /// \f\[
    ///     b'(d) =
    ///     2\frac{1}{\hat{d}}\left(1-\frac{d}{\hat{d}}\right)\ln\left(\frac{d}{\hat{d}}\right)
    ///             + \left(1-\frac{d}{\hat{d}}\right)^2 \frac{1}{d}
    /// \f\]
    ///
    /// @param d The distance.
    /// @param dhat Activation distance of the barrier.
    /// @return The derivative of the barrier wrt d.
    T first_derivative(const T d, const T dhat) const override
    {
        return BarrierT::first_derivative(d / dhat, T(1)) / dhat;
    }

    /// @brief Second derivative of the barrier function.
    ///
    /// \f\[
    ///     b''(d) = \frac{\hat{d}^2-2 d^2 \ln \left(\frac{d}{\hat{d}}\right)+2
    ///     \hat{d} d-3 d^2}{\hat{d}^2 d^2}
    /// \f\]
    ///
    /// @param d The distance.
    /// @param dhat Activation distance of the barrier.
    /// @return The second derivative of the barrier wrt d.
    T second_derivative(const T d, const T dhat) const override
    {
        return BarrierT::second_derivative(d / dhat, T(1)) / (dhat * dhat);
    }

    /// @brief Get the units of the barrier function.
    /// @param dhat The activation distance of the barrier.
    /// @return The units of the barrier function.
    T units(const T dhat) const override
    {
        return T(1); // The normalized barrier is dimensionless.
    }
};

template <typename T = double>
using NormalizedClampedLogBarrier = NormalizedBarrier<ClampedLogBarrier<T>>;

// ============================================================================
// Quadratic log barrier functions from [Huang et al. 2024]
// ============================================================================

/// @brief Clamped log barrier with a quadratic log term from [Huang et al. 2024].
template <typename T = double>
class ClampedLogSqBarrier : public BarrierBase<T> {
public:
    ClampedLogSqBarrier() = default;

    /// @brief Function that grows to infinity as d approaches 0 from the right.
    ///
    /// \f\[
    ///     b(d) = (d-\hat{d})^2\ln^2\left(\frac{d}{\hat{d}}\right)
    /// \f\]
    ///
    /// @param d The distance.
    /// @param dhat Activation distance of the barrier.
    /// @return The value of the barrier function at d.
    T operator()(const T d, const T dhat) const override;

    /// @brief Derivative of the barrier function.
    ///
    /// \f\[
    ///     b'(d) = 2 (d - \hat{d}) \ln\left(\frac{d}{\hat{d}}\right)
    ///     \left[\ln\left(\frac{d}{\hat{d}}\right) + \frac{d -
    ///     \hat{d}}{d}\right]
    /// \f\]
    ///
    /// @param d The distance.
    /// @param dhat Activation distance of the barrier.
    /// @return The derivative of the barrier wrt d.
    T first_derivative(const T d, const T dhat) const override;

    /// @brief Second derivative of the barrier function.
    ///
    /// \f\[
    ///     b''(d) = 2 \left(\ln^2\left(\frac{d}{\hat{d}}\right) - \left(
    ///     \ln\left(\frac{d}{\hat{d}}\right) - 1\right) \frac{(\hat{d} -
    ///     d)^2}{d^2} - 4 \ln\left(\frac{d}{\hat{d}}\right) \frac{\hat{d} -
    ///     d}{d}\right)
    /// \f\]
    ///
    /// @param d The distance.
    /// @param dhat Activation distance of the barrier.
    /// @return The second derivative of the barrier wrt d.
    T second_derivative(const T d, const T dhat) const override;

    /// @brief Get the units of the barrier function.
    /// @param dhat The activation distance of the barrier.
    /// @return The units of the barrier function.
    T units(const T dhat) const override
    {
        // (d - d̂)² = d̂² (d/d̂ - 1)²
        return dhat * dhat;
    }
};

// ============================================================================
// Cubic barrier from [Ando 2024]
// ============================================================================

/// @brief Cubic barrier function from [Ando 2024].
template <typename T = double> class CubicBarrier : public BarrierBase<T> {
public:
    CubicBarrier() = default;

    /// @brief Weak barrier function.
    ///
    /// \f\[
    ///     b(d) = -\frac{2}{3\hat{d}} (d - \hat{d})^3
    /// \f\]
    ///
    /// @param d The distance.
    /// @param dhat Activation distance of the barrier.
    /// @return The value of the barrier function at d.
    T operator()(const T d, const T dhat) const override;

    /// @brief Derivative of the barrier function.
    ///
    /// \f\[
    ///     b'(d) = -2 (d - \hat{d})^2
    /// \f\]
    ///
    /// @param d The distance.
    /// @param dhat Activation distance of the barrier.
    /// @return The derivative of the barrier wrt d.
    T first_derivative(const T d, const T dhat) const override;

    /// @brief Second derivative of the barrier function.
    ///
    /// \f\[
    ///     b''(d) = -4 (d - \hat{d})
    /// \f\]
    ///
    /// @param d The distance.
    /// @param dhat Activation distance of the barrier.
    /// @return The second derivative of the barrier wrt d.
    T second_derivative(const T d, const T dhat) const override;

    /// @brief Get the units of the barrier function.
    /// @param dhat The activation distance of the barrier.
    /// @return The units of the barrier function.
    T units(const T dhat) const override
    {
        // (d - d̂)² = d̂² (d/d̂ - 1)²
        return dhat * dhat;
    }
};

// ============================================================================
// 2-Stage activation function from [Chen et al. 2025]
// ============================================================================

/// @brief 2-Stage activation function from [Chen et al. 2025].
template <typename T = double> class TwoStageBarrier : public BarrierBase<T> {
public:
    TwoStageBarrier() = default;

    /**
     * @brief Two-stage activation barrier.
     *
     * \f\[
     *     b(d) = \begin{cases}
     *         -\frac{\hat{d}^2}{4} \left(\ln\left(\frac{2d}{\hat{d}}\right) -
     *         \tfrac{1}{2}\right) & d < \frac{\hat{d}}{2}\\
     *         \tfrac{1}{2} (\hat{d} - d)^2 & d < \hat{d}\\
     *         0 & d \ge \hat{d}
     *     \end{cases}
     * \f\]
     *
     * @param d The distance.
     * @param dhat Activation distance of the barrier.
     * @return The value of the barrier function at d.
     */
    T operator()(const T d, const T dhat) const override;

    /**
     * @brief Derivative of the barrier function.
     *
     * \f\[
     *     b'(d) = \begin{cases}
     *         -\frac{\hat{d}}{4d} & d < \frac{\hat{d}}{2}\\
     *         d - \hat{d} & d < \hat{d}\\
     *         0 & d \ge \hat{d}
     *     \end{cases}
     * \f\]
     *
     * @param d The distance.
     * @param dhat Activation distance of the barrier.
     * @return The derivative of the barrier wrt d.
     */
    T first_derivative(const T d, const T dhat) const override;

    /**
     * @brief Second derivative of the barrier function.
     *
     * \f\[
     *     b''(d) = \begin{cases}
     *         \frac{\hat{d}}{4d^2} & d < \frac{\hat{d}}{2}\\
     *         1 & d < \hat{d}\\
     *         0 & d \ge \hat{d}
     *     \end{cases}
     * \f\]
     *
     * @param d The distance.
     * @param dhat Activation distance of the barrier.
     * @return The second derivative of the barrier wrt d.
     */
    T second_derivative(const T d, const T dhat) const override;

    /// @brief Get the units of the barrier function.
    /// @param dhat The activation distance of the barrier.
    /// @return The units of the barrier function.
    T units(const T dhat) const override
    {
        // (d - d̂)² = d̂² (d/d̂ - 1)²
        return dhat * dhat;
    }
};

} // namespace ipc
