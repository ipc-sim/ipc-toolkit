// Barrier functions that grow to infinity as x -> 0+. Includes gradient and
// hessian functions, too. These barrier functions can be used to impose
// inequality constraints on a function.

#pragma once

#include <ipc/config.hpp>

#include <cmath>
#include <memory>

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
template <typename T = double>
IPC_TOOLKIT_HOST_DEVICE T barrier(const T d, const T dhat);

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
IPC_TOOLKIT_HOST_DEVICE T barrier_first_derivative(const T d, const T dhat);

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
IPC_TOOLKIT_HOST_DEVICE T barrier_second_derivative(const T d, const T dhat);

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
    using Scalar = typename BarrierT::value_type;

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
    Scalar operator()(const Scalar d, const Scalar dhat) const override
    {
        return BarrierT::operator()(d / dhat, Scalar(1));
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
    Scalar first_derivative(const Scalar d, const Scalar dhat) const override
    {
        return BarrierT::first_derivative(d / dhat, Scalar(1)) / dhat;
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
    Scalar second_derivative(const Scalar d, const Scalar dhat) const override
    {
        return BarrierT::second_derivative(d / dhat, Scalar(1)) / (dhat * dhat);
    }

    /// @brief Get the units of the barrier function.
    /// @param dhat The activation distance of the barrier.
    /// @return The units of the barrier function.
    Scalar units(const Scalar dhat) const override
    {
        return Scalar(1); // The normalized barrier is dimensionless.
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
     *         \infty & d \le 0\\
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
     *         0 & d \le 0\\
     *         -\frac{\hat{d}^2}{4d} & d < \frac{\hat{d}}{2}\\
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
     *         0 & d \le 0\\
     *         \frac{\hat{d}^2}{4d^2} & d < \frac{\hat{d}}{2}\\
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

// ============================================================================
// Inverse-power barrier
// ============================================================================

/// @brief Inverse-power barrier with smooth compact support.
///
/// \f\[
///     b(d) = \frac{h(d,\hat{d})}{d^p}, \quad
///     h(d,\hat{d}) = 2\,B\!\left(\frac{2d}{\hat{d}}\right)
/// \f\]
///
/// where \f$B\f$ is the standard cubic B-spline basis function and \f$p > 0\f$
/// is the power parameter. The window \f$h\f$ vanishes smoothly at
/// \f$d = \hat{d}\f$ (C² contact), ensuring \f$b(d)=0\f$ for \f$d\ge\hat{d}\f$,
/// while \f$b(d)\to+\infty\f$ as \f$d\to 0^+\f$.
class InversePowerBarrier : public Barrier {
public:
    /// @param power The power \f$p > 0\f$ controlling the singularity at d = 0.
    explicit InversePowerBarrier(const double power) : m_power(power) { }

    /// @brief b(d) = h(d, d̂) / d^p.
    /// @param d Distance (must be > 0 for a finite value).
    /// @param dhat Activation distance of the barrier.
    /// @return The value of the barrier function at d.
    double operator()(const double d, const double dhat) const override;

    /// @brief First derivative b'(d) = (h'·d − p·h) / d^(p+1).
    /// @param d Distance.
    /// @param dhat Activation distance of the barrier.
    /// @return The first derivative of the barrier function at d.
    double first_derivative(const double d, const double dhat) const override;

    /// @brief Second derivative b''(d) = (h''·d² − 2p·h'·d + p(p+1)·h) / d^(p+2).
    /// @param d Distance.
    /// @param dhat Activation distance of the barrier.
    /// @return The second derivative of the barrier function at d.
    double second_derivative(const double d, const double dhat) const override;

    /// @brief Get the units of the barrier function (d̂^{-p}).
    /// @param dhat The activation distance of the barrier.
    /// @return The units of the barrier function.
    double units(const double dhat) const override
    {
        return 1.0 / std::pow(dhat, m_power);
    }

    /// @brief The power p used in the barrier.
    double power() const { return m_power; }

private:
    double m_power; ///< p > 0

    /// @brief Evaluate the B-spline window h(d, dhat) and its first two
    ///        derivatives with respect to d.
    ///
    /// h(d) = 2 * B(2d/dhat) where B is the cubic B-spline:
    ///   B(t) = 2/3 - t² + t³/2        for 0 ≤ t < 1
    ///   B(t) = (2-t)³ / 6             for 1 ≤ t < 2
    ///   B(t) = 0                       for t ≥ 2
    static void
    h_and_derivs(double d, double dhat, double& h, double& dh, double& ddh);
};

/// @brief Near-Far barrier function.
/// This barrier function takes another "base" barrier object as argument.
/// All operations defined in Barrier defer the computation to its base barrier.
class NearFarBarrier : public Barrier {
public:
    /// @brief Construct a NearFarBarrier.
    /// @param base_barrier The base barrier function to use.
    /// @param alpha A double parameter.
    NearFarBarrier(const Barrier* const base_barrier, const double alpha)
        : m_base_barrier(base_barrier)
        , m_alpha(alpha)
    {
    }

    /// @brief Construct a NearFarBarrier holding shared ownership of the base
    /// barrier. Prefer this over the raw-pointer overload whenever the caller
    /// already has a shared_ptr — it pins the base barrier for the lifetime of
    /// this object, eliminating dangling-pointer risk.
    NearFarBarrier(
        std::shared_ptr<const Barrier> base_barrier, const double alpha)
        : m_base_barrier_owned(std::move(base_barrier))
        , m_base_barrier(m_base_barrier_owned.get())
        , m_alpha(alpha)
    {
    }

    /// @brief Evaluate the barrier function.
    /// @param d Distance.
    /// @param dhat Activation distance of the barrier.
    /// @return The value of the barrier function at d.
    double operator()(const double d, const double dhat) const override
    {
        return (*m_base_barrier)(d, dhat);
    }

    /// @brief Evaluate the first derivative of the barrier function wrt d.
    /// @param d Distance.
    /// @param dhat Activation distance of the barrier.
    /// @return The value of the first derivative of the barrier function at d.
    double first_derivative(const double d, const double dhat) const override
    {
        return m_base_barrier->first_derivative(d, dhat);
    }

    /// @brief Evaluate the second derivative of the barrier function wrt d.
    /// @param d Distance.
    /// @param dhat Activation distance of the barrier.
    /// @return The value of the second derivative of the barrier function at d.
    double second_derivative(const double d, const double dhat) const override
    {
        return m_base_barrier->second_derivative(d, dhat);
    }

    /// @brief Get the units of the barrier function.
    /// @param dhat The activation distance of the barrier.
    /// @return The units of the barrier function.
    double units(const double dhat) const override
    {
        return m_base_barrier->units(dhat);
    }

    /// @brief Evaluate the near function.
    double near_value(const double d, const double dhat) const;

    /// @brief Evaluate the far function.
    double far_value(const double d, const double dhat) const;

    /// @brief Evaluate the first derivative of the near function.
    double first_derivative_near(const double d, const double dhat) const;

    /// @brief Evaluate the first derivative of the far function.
    double first_derivative_far(const double d, const double dhat) const;

    /// @brief Evaluate the second derivative of the near function.
    double second_derivative_near(const double d, const double dhat) const;

    /// @brief Evaluate the second derivative of the far function.
    double second_derivative_far(const double d, const double dhat) const;

private:
    // Optional shared ownership; null when constructed from a raw pointer.
    const std::shared_ptr<const Barrier> m_base_barrier_owned;
    const Barrier* const m_base_barrier;
    const double m_alpha;
};

} // namespace ipc
