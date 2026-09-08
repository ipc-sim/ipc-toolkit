// Out-of-line definitions for the Barrier class hierarchy.
//
// These live apart from barrier.cpp because that file is compiled as a CUDA
// translation unit when CUDA is enabled, and a virtual hierarchy has nothing to
// offer device code anyway. Keeping the classes here lets barrier.cpp stay
// small and device-callable, holding only the three free barrier functions.
#include "barrier.hpp"

#include <ipc/config.hpp>
#include <ipc/math/scalar_math.hpp>
#include <ipc/utils/simd.hpp>

namespace ipc {

// ============================================================================

template <typename T>
T ClampedLogSqBarrier<T>::operator()(const T d, const T dhat) const
{
    using namespace ipc::numext; // log
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
    using namespace ipc::numext; // log
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
    using namespace ipc::numext; // log
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
    using namespace ipc::numext; // log
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
#endif
/// @endcond
// ============================================================================

} // namespace ipc
