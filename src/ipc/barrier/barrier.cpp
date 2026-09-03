#include <ipc/config.hpp>

// Barrier functions that grow to infinity as x -> 0+. Includes gradient and
// hessian functions, too. These barrier functions can be used to impose
// inequality constraints on a function.
#include "barrier.hpp"

#include <cmath>
#include <limits>

namespace ipc {

template <typename T> T barrier(const T d, const T dhat)
{
    if (d <= T(0)) {
        return std::numeric_limits<T>::infinity();
    }
    if (d >= dhat) {
        return T(0);
    }
    // b(d) = -(d-d̂)²ln(d / d̂)
    const T d_minus_dhat = (d - dhat);
    return -d_minus_dhat * d_minus_dhat * std::log(d / dhat);
}

template <typename T> T barrier_first_derivative(const T d, const T dhat)
{
    if (d <= T(0) || d >= dhat) {
        return T(0);
    }
    // b(d) = -(d - d̂)²ln(d / d̂)
    // b'(d) = -2(d - d̂)ln(d / d̂) - (d-d̂)²(1 / d)
    //       = (d - d̂) * (-2ln(d/d̂) - (d - d̂) / d)
    //       = (d̂ - d) * (2ln(d/d̂) - d̂/d + 1)
    return (dhat - d) * (2 * std::log(d / dhat) - dhat / d + 1);
}

template <typename T> T barrier_second_derivative(const T d, const T dhat)
{
    if (d <= T(0) || d >= dhat) {
        return T(0);
    }
    const T dhat_d = dhat / d;
    return (dhat_d + 2) * dhat_d - 2 * std::log(d / dhat) - 3;
}

// ============================================================================

template <typename T>
T ClampedLogSqBarrier<T>::operator()(const T d, const T dhat) const
{
    if (d <= T(0)) {
        return std::numeric_limits<T>::infinity();
    }
    if (d >= dhat) {
        return T(0);
    }
    // b(d) = (d-d̂)²ln²(d / d̂)
    const T d_minus_dhat = (d - dhat);
    const T log_d_dhat = std::log(d / dhat);
    return d_minus_dhat * d_minus_dhat * log_d_dhat * log_d_dhat;
}

template <typename T>
T ClampedLogSqBarrier<T>::first_derivative(const T d, const T dhat) const
{
    if (d <= T(0) || d >= dhat) {
        return T(0);
    }
    // b(d) = (d - d̂)²ln²(d / d̂)
    // b'(d) = 2 (d - d̂) ln²(d / d̂) + 2 (d - d̂)² ln(d / d̂) / d
    //       = 2 (d - d̂) ln(d / d̂) [ln(d / d̂) + (d - d̂) / d]
    const T d_minus_dhat = (d - dhat);
    const T log_d_dhat = std::log(d / dhat);
    return T(2) * d_minus_dhat * log_d_dhat * (log_d_dhat + d_minus_dhat / d);
}

template <typename T>
T ClampedLogSqBarrier<T>::second_derivative(const T d, const T dhat) const
{
    if (d <= T(0) || d >= dhat) {
        return T(0);
    }
    const T t0 = dhat - d;
    const T t1 = std::log(d / dhat);
    const T t2 = (t0 * t0) / (d * d);
    return T(2) * ((t1 * t1) - (t1 - T(1)) * t2 - T(4) * t1 * t0 / d);
}

// ============================================================================

template <typename T>
T CubicBarrier<T>::operator()(const T d, const T dhat) const
{
    if (d < dhat) {
        // b(d) = (d - d̂)³
        const T d_minus_dhat = (d - dhat);
        return -T(2) / T(3) / dhat * d_minus_dhat * d_minus_dhat * d_minus_dhat;
    } else {
        return T(0);
    }
}

template <typename T>
T CubicBarrier<T>::first_derivative(const T d, const T dhat) const
{
    if (d < dhat) {
        const T d_minus_dhat = (d - dhat);
        return T(-2) / dhat * d_minus_dhat * d_minus_dhat;
    } else {
        return T(0);
    }
}

template <typename T>
T CubicBarrier<T>::second_derivative(const T d, const T dhat) const
{
    if (d < dhat) {
        return T(4) * (T(1) - d / dhat);
    } else {
        return T(0);
    }
}

// ============================================================================

template <typename T>
T TwoStageBarrier<T>::operator()(const T d, const T dhat) const
{
    if (d >= dhat) {
        return T(0);
    } else if (d >= T(0.5) * dhat) {
        return T(0.5) * (dhat - d) * (dhat - d);
    } else {
        return T(-0.25) * dhat * dhat * (std::log(T(2) * d / dhat) - T(0.5));
    }
}

template <typename T>
T TwoStageBarrier<T>::first_derivative(const T d, const T dhat) const
{
    if (d >= dhat) {
        return T(0);
    } else if (d >= T(0.5) * dhat) {
        return d - dhat;
    } else {
        return T(-0.25) * dhat * dhat / d;
    }
}

template <typename T>
T TwoStageBarrier<T>::second_derivative(const T d, const T dhat) const
{
    if (d >= dhat) {
        return T(0);
    } else if (d >= T(0.5) * dhat) {
        return T(1);
    } else {
        return (T(0.25) * dhat * dhat) / (d * d);
    }
}

// ============================================================================
// Explicit template instantiations
/// @cond DOXYGEN_SKIP
template class BarrierBase<float>;
template class BarrierBase<double>;
template class ClampedLogBarrier<float>;
template class ClampedLogBarrier<double>;
template class ClampedLogSqBarrier<float>;
template class ClampedLogSqBarrier<double>;
template class CubicBarrier<float>;
template class CubicBarrier<double>;
template class TwoStageBarrier<float>;
template class TwoStageBarrier<double>;
template float barrier(const float d, const float dhat);
template double barrier(const double d, const double dhat);
template float barrier_first_derivative(const float d, const float dhat);
template double barrier_first_derivative(const double d, const double dhat);
template float barrier_second_derivative(const float d, const float dhat);
template double barrier_second_derivative(const double d, const double dhat);
/// @endcond
// ============================================================================

} // namespace ipc
