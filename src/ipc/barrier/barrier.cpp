#include <ipc/config.hpp>

// Barrier functions that grow to infinity as x -> 0+. Includes gradient and
// hessian functions, too. These barrier functions can be used to impose
// inequality constraints on a function.
#include "barrier.hpp"

#include <ipc/math/scalar_math.hpp>
#include <ipc/utils/simd.hpp>

namespace ipc {

// Each barrier is one select_lazy cascade, ordered by increasing d so it
// reads like the piecewise definition in the header. A scalar evaluates only
// the case it lands in -- so the log below is never reached for d <= 0 -- while
// a batch evaluates every case and blends per-lane, earlier cases winning.

template <typename T> T barrier(const T d, const T dhat)
{
    // b(d) = -(d-d̂)²ln(d / d̂)
    return select_lazy(
        d <= T(0), [&] { return infinity<T>(); }, //
        d < dhat, [&] { return -sqr(d - dhat) * log(d / dhat); },
        [&] { return T(0); });
}

template <typename T> T barrier_first_derivative(const T d, const T dhat)
{
    // b(d) = -(d - d̂)²ln(d / d̂)
    // b'(d) = -2(d - d̂)ln(d / d̂) - (d-d̂)²(1 / d)
    //       = (d - d̂) * (-2ln(d/d̂) - (d - d̂) / d)
    //       = (d̂ - d) * (2ln(d/d̂) - d̂/d + 1)
    return select_lazy(
        d <= T(0), [&] { return T(0); }, //
        d < dhat,
        [&] { return (dhat - d) * (2 * log(d / dhat) - dhat / d + 1); },
        [&] { return T(0); });
}

template <typename T> T barrier_second_derivative(const T d, const T dhat)
{
    return select_lazy(
        d <= T(0), [&] { return T(0); }, //
        d < dhat,
        [&] {
            const T dhat_d = dhat / d;
            return (dhat_d + 2) * dhat_d - 2 * log(d / dhat) - 3;
        },
        [&] { return T(0); });
}

// ============================================================================

template <typename T>
T ClampedLogSqBarrier<T>::operator()(const T d, const T dhat) const
{
    // b(d) = (d-d̂)²ln²(d / d̂)
    return select_lazy(
        d <= T(0), [&] { return infinity<T>(); }, //
        d < dhat,
        [&] {
            const T log_d_dhat = log(d / dhat);
            return sqr(d - dhat) * sqr(log_d_dhat);
        },
        [&] { return T(0); });
}

template <typename T>
T ClampedLogSqBarrier<T>::first_derivative(const T d, const T dhat) const
{
    // b(d) = (d - d̂)²ln²(d / d̂)
    // b'(d) = 2 (d - d̂) ln²(d / d̂) + 2 (d - d̂)² ln(d / d̂) / d
    //       = 2 (d - d̂) ln(d / d̂) [ln(d / d̂) + (d - d̂) / d]
    return select_lazy(
        d <= T(0), [&] { return T(0); }, //
        d < dhat,
        [&] {
            const T d_minus_dhat = (d - dhat);
            const T log_d_dhat = log(d / dhat);
            return T(2) * d_minus_dhat * log_d_dhat
                * (log_d_dhat + d_minus_dhat / d);
        },
        [&] { return T(0); });
}

template <typename T>
T ClampedLogSqBarrier<T>::second_derivative(const T d, const T dhat) const
{
    return select_lazy(
        d <= T(0), [&] { return T(0); }, //
        d < dhat,
        [&] {
            const T t0 = dhat - d;
            const T t1 = log(d / dhat);
            const T t2 = sqr(t0) / sqr(d);
            return T(2) * (sqr(t1) - (t1 - T(1)) * t2 - T(4) * t1 * t0 / d);
        },
        [&] { return T(0); });
}

// ============================================================================

template <typename T>
T CubicBarrier<T>::operator()(const T d, const T dhat) const
{
    // b(d) = (d - d̂)³
    //
    // A polynomial: finite at d <= 0, so unlike the log barriers it needs no
    // penetration guard.
    return select_lazy(
        d < dhat, [&] { return -T(2) / T(3) / dhat * cubic(d - dhat); },
        [&] { return T(0); });
}

template <typename T>
T CubicBarrier<T>::first_derivative(const T d, const T dhat) const
{
    return select_lazy(
        d < dhat, [&] { return T(-2) / dhat * sqr(d - dhat); },
        [&] { return T(0); });
}

template <typename T>
T CubicBarrier<T>::second_derivative(const T d, const T dhat) const
{
    return select_lazy(
        d < dhat, [&] { return T(4) * (T(1) - d / dhat); },
        [&] { return T(0); });
}

// ============================================================================

template <typename T>
T TwoStageBarrier<T>::operator()(const T d, const T dhat) const
{
    return select_lazy(
        d <= T(0), [&] { return infinity<T>(); }, //
        d < T(0.5) * dhat,
        [&] { return T(-0.25) * sqr(dhat) * (log(T(2) * d / dhat) - T(0.5)); },
        d < dhat, [&] { return T(0.5) * sqr(dhat - d); }, //
        [&] { return T(0); });
}

template <typename T>
T TwoStageBarrier<T>::first_derivative(const T d, const T dhat) const
{
    return select_lazy(
        d <= T(0), [&] { return T(0); },                             //
        d < T(0.5) * dhat, [&] { return T(-0.25) * sqr(dhat) / d; }, //
        d < dhat, [&] { return d - dhat; },                          //
        [&] { return T(0); });
}

template <typename T>
T TwoStageBarrier<T>::second_derivative(const T d, const T dhat) const
{
    return select_lazy(
        d <= T(0), [&] { return T(0); },                                   //
        d < T(0.5) * dhat, [&] { return (T(0.25) * sqr(dhat)) / sqr(d); }, //
        d < dhat, [&] { return T(1); },                                    //
        [&] { return T(0); });
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
#ifdef IPC_TOOLKIT_WITH_SIMD
template class BarrierBase<SimdBatch<float>>;
template class BarrierBase<SimdBatch<double>>;
template class ClampedLogBarrier<SimdBatch<float>>;
template class ClampedLogBarrier<SimdBatch<double>>;
template class ClampedLogSqBarrier<SimdBatch<float>>;
template class ClampedLogSqBarrier<SimdBatch<double>>;
template class CubicBarrier<SimdBatch<float>>;
template class CubicBarrier<SimdBatch<double>>;
template class TwoStageBarrier<SimdBatch<float>>;
template class TwoStageBarrier<SimdBatch<double>>;
template SimdBatch<float>
barrier(const SimdBatch<float> d, const SimdBatch<float> dhat);
template SimdBatch<double>
barrier(const SimdBatch<double> d, const SimdBatch<double> dhat);
template SimdBatch<float>
barrier_first_derivative(const SimdBatch<float> d, const SimdBatch<float> dhat);
template SimdBatch<double> barrier_first_derivative(
    const SimdBatch<double> d, const SimdBatch<double> dhat);
template SimdBatch<float> barrier_second_derivative(
    const SimdBatch<float> d, const SimdBatch<float> dhat);
template SimdBatch<double> barrier_second_derivative(
    const SimdBatch<double> d, const SimdBatch<double> dhat);
#endif
/// @endcond
// ============================================================================

} // namespace ipc
