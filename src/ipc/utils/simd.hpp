#pragma once

#include <ipc/config.hpp>
#include <ipc/math/scalar_math.hpp>

#include <Eigen/Core>

#include <limits>
#include <type_traits>
#include <utility>

namespace ipc {

/// @brief Whether `T` is an `xsimd::batch`, i.e. a per-lane scalar with no
/// single `bool` value: comparisons and reductions like `isZero()` cannot be
/// used in scalar control flow (`if`, `assert`) for such a `T`.
///
/// Specialized below, where the batch type itself is available.
template <typename T> inline constexpr bool is_simd_batch_v = false;

/// @brief Pick between `a` and `b`.
///
/// The scalar counterpart of `xsimd::select`, which ADL finds for a batch
/// `mask`, so one `select(cond, a, b)` compiles for both.
template <typename T> inline T select(const bool mask, const T& a, const T& b)
{
    return mask ? a : b;
}

/// @brief `+infinity` for any scalar the library templates on.
template <typename T> inline T infinity()
{
    if constexpr (is_simd_batch_v<T>) {
        return T(std::numeric_limits<typename T::value_type>::infinity());
    } else {
        return std::numeric_limits<T>::infinity();
    }
}

/// @brief A first-match-wins cascade of cases:
/// `select_lazy(mask1, value1, mask2, value2, ..., else_value)`.
///
/// Each value is a callable, which is what lets the two scalar/batch
/// semantics share one definition:
///
/// - For a scalar `T`, this is an `if`/`else if` chain: only the first
///   matching value is evaluated. A case guarded by `d > 0` can therefore
///   hold a `log(d)` safely.
/// - For a batch, no mask has one answer, so *every* value is evaluated and
///   blended per-lane. One batch may carry lanes in different cases. A value
///   may be evaluated on lanes its mask excludes -- producing a NaN that is
///   then blended away -- so it must not have side effects.
///
/// The masks are ordinary arguments, so unlike the values they are all
/// evaluated before the call, whichever case wins. Keep them to comparisons.
///
/// @warning Earlier cases take precedence in both paths, which for the batch
/// comes from the order the blend is folded. Masks may overlap, and only that
/// order decides the winner, so a test must cover lanes that fall in
/// overlapping cases.
template <typename F> inline auto select_lazy(F&& else_value)
{
    return else_value();
}

template <typename Mask, typename F, typename... Rest>
inline auto select_lazy(const Mask& mask, F&& value, Rest&&... rest)
{
    if constexpr (std::is_same_v<std::decay_t<Mask>, bool>) {
        return mask ? value() : select_lazy(std::forward<Rest>(rest)...);
    } else {
        return select(mask, value(), select_lazy(std::forward<Rest>(rest)...));
    }
}

} // namespace ipc

#ifdef IPC_TOOLKIT_WITH_SIMD

#include <xsimd/xsimd.hpp>

namespace Eigen {

/// @brief Let Eigen treat an xsimd batch as a scalar type.
///
/// This makes ``Eigen::Vector3<ipc::SimdBatch<double>>`` a valid argument to
/// the distance functions, evaluating one independent problem per SIMD lane.
template <typename T, typename A>
struct NumTraits<xsimd::batch<T, A>> : GenericNumTraits<xsimd::batch<T, A>> {
    using Real = xsimd::batch<T, A>;
    using NonInteger = xsimd::batch<T, A>;
    using Nested = xsimd::batch<T, A>;
    using Literal = xsimd::batch<T, A>;

    // NOLINTBEGIN(readability-identifier-naming,performance-enum-size)
    enum {
        IsComplex = 0,
        IsInteger = 0,
        IsSigned = 1,
        // A batch is not trivially default-constructible in the sense Eigen
        // means; make it initialize its coefficients.
        RequireInitialization = 1,
        ReadCost = 1,
        AddCost = 1,
        MulCost = 1
    };
    // NOLINTEND(readability-identifier-naming,performance-enum-size)

    static Real epsilon() { return Real(NumTraits<T>::epsilon()); }
    static Real dummy_precision()
    {
        return Real(NumTraits<T>::dummy_precision());
    }
    static int digits10() { return NumTraits<T>::digits10(); }
    static Real highest() { return Real(NumTraits<T>::highest()); }
    static Real lowest() { return Real(NumTraits<T>::lowest()); }
};

} // namespace Eigen

namespace ipc {

/// @brief A SIMD batch of ``T``, usable as the scalar type of a distance query.
///
/// Packing one independent problem per lane (a structure-of-arrays layout)
/// lets a single call evaluate ``SimdBatch<T>::size`` problems at once:
///
/// @code
/// Eigen::Vector3<ipc::SimdBatch<double>> p, e0, e1;
/// for (int k = 0; k < 3; ++k) {
///     p[k] = ipc::SimdBatch<double>::load_unaligned(&p_soa[k * n + i]);
///     // ... likewise for e0, e1
/// }
/// const auto d = ipc::point_edge_distance(p, e0, e1,
/// PointEdgeDistanceType::P_E);
/// @endcode
///
/// @warning **An explicit distance type is required.** ``AUTO`` and the
/// ``*_distance_type`` predicates are not available for batch scalars, and
/// throw ``std::invalid_argument`` if reached: the distance type is a per-lane
/// property, but the predicates return a single enum, so two lanes of the same
/// batch cannot report different closest features. Resolve the distance types
/// scalar-side and group the problems by type before batching, which is the
/// natural structure-of-arrays layout anyway.
///
/// @note Functions that select a case from the values themselves — clamping a
/// barrier at ``d̂``, choosing a tangent-basis reference axis — evaluate every
/// case for a batch and blend the results per-lane, so a batch can carry lanes
/// that fall in different cases. The ``*_distance_type`` predicates are the
/// exception above: their result is an enum, not a number to blend.
///
/// @warning ``xsimd::default_arch`` is chosen from the compiler flags of each
/// translation unit, so a caller compiled without the SIMD flags the library
/// was built with names a *different* type than the one instantiated here, and
/// will fail to link. Build such callers with the same flags — CMake reports
/// them as ``SIMD_CXX_FLAGS``.
template <typename T> using SimdBatch = xsimd::batch<T, xsimd::default_arch>;

template <typename T, typename A>
inline constexpr bool is_simd_batch_v<xsimd::batch<T, A>> = true;

} // namespace ipc

#endif

namespace ipc {

/// @brief `v.normalized()`, but also defined for a batch scalar.
///
/// Eigen's own `normalized()` leaves a zero-length vector unscaled behind an
/// `if (squaredNorm() > 0)`, which a batch cannot answer with one bool. This
/// applies that same rule per-lane and otherwise defers to Eigen.
template <typename Derived>
inline typename Derived::PlainObject
normalized(const Eigen::MatrixBase<Derived>& v)
{
    using T = typename Derived::Scalar;
    if constexpr (is_simd_batch_v<T>) {
        const typename Derived::PlainObject n = v.derived();
        const T z = n.squaredNorm();
        const typename Derived::PlainObject scaled = n / ipc::sqrt(z);
        const auto is_nonzero = z > T(0);

        typename Derived::PlainObject out = n;
        for (Eigen::Index i = 0; i < out.size(); ++i) {
            out(i) = select(is_nonzero, scaled(i), n(i));
        }
        return out;
    } else {
        return v.normalized();
    }
}

} // namespace ipc
