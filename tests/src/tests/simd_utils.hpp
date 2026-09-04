#pragma once

#include <ipc/utils/simd.hpp>

#ifdef IPC_TOOLKIT_WITH_SIMD

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <Eigen/Core>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <string>

namespace ipc::tests {

/// @brief The batch type the SIMD tests pack their lanes into.
using Batch = ipc::SimdBatch<double>;

/// @brief Lanes per batch: the number of independent problems one call solves.
constexpr int L = int(Batch::size);

/// @brief The alignment a batch load/store wants.
///
/// The scalar fallback architecture reports 0, which `alignas` treats as "no
/// effect" but which would under-align the array on its own, so take the
/// scalar's own alignment as the floor.
constexpr std::size_t LANE_ALIGNMENT = Batch::arch_type::alignment()
        > alignof(double)
    ? Batch::arch_type::alignment()
    : alignof(double);

/// @brief One lane's worth of scalars per lane, aligned for a batch
/// load/store.
///
/// The `alignas` is what makes `load_aligned`/`store_aligned` usable: a plain
/// `std::array<double, L>` is only guaranteed `alignof(double)`, and an
/// aligned load off an under-aligned address faults on the architectures that
/// require alignment rather than merely running slower.
struct alignas(LANE_ALIGNMENT) Lanes {
    std::array<double, L> values {};

    double& operator[](const int l) { return values[size_t(l)]; }
    const double& operator[](const int l) const { return values[size_t(l)]; }
    double* data() { return values.data(); }
    const double* data() const { return values.data(); }

    /// @brief Load these lanes into one batch.
    Batch load() const { return Batch::load_aligned(data()); }

    /// @brief Unpack a batch back into its lanes.
    static Lanes from(const Batch& b)
    {
        Lanes out;
        b.store_aligned(out.data());
        return out;
    }
};

static_assert(
    alignof(Lanes) >= LANE_ALIGNMENT,
    "Lanes must satisfy the alignment its load/store claims");

/// @brief `L` problems, one per lane.
template <int dim> using Points = std::array<Eigen::Vector<double, dim>, L>;

// ---------------------------------------------------------------------------
// Tolerances

/// @brief Agreement bound for a value.
///
/// The batch and scalar paths differ only in how the compiler contracts
/// multiply-adds, which for a value costs no more than a rounding step.
constexpr double VALUE_TOL = 1e-14;

/// @brief Agreement bound for a derivative.
///
/// Applied relative to the magnitude of the *whole* result rather than
/// per-entry: a gradient or Hessian of a near-degenerate configuration has
/// entries that are catastrophic cancellations of much larger terms, so a
/// per-entry bound would reject a correct result. A structurally wrong entry
/// differs by order of the result's own magnitude.
constexpr double DERIVATIVE_TOL = 1e-9;

// ---------------------------------------------------------------------------
// Packing

/// @brief Pack the k-th coordinate of `L` problems into one batch each.
template <int dim> inline Eigen::Vector<Batch, dim> pack(const Points<dim>& v)
{
    Eigen::Vector<Batch, dim> out;
    for (int k = 0; k < dim; ++k) {
        Lanes lane;
        for (int l = 0; l < L; ++l) {
            lane[l] = v[size_t(l)][k];
        }
        out[k] = lane.load();
    }
    return out;
}

/// @brief `L` random points, deterministic in `seed`.
template <int dim> inline Points<dim> random_points(const int seed)
{
    std::srand(unsigned(seed));
    Points<dim> v;
    for (int l = 0; l < L; ++l) {
        v[size_t(l)] = Eigen::Vector<double, dim>::Random();
    }
    return v;
}

/// @brief Assign `cases` round-robin to the `L` lanes, starting at `offset`.
///
/// A batch may hold fewer lanes than there are cases, so a test sweeps the
/// offset over `cases.size()` to reach every case — and, where the two counts
/// differ, to land different cases in the same batch, which is what exercises
/// the per-lane blend.
template <typename Container>
inline Lanes lane_cases(const Container& cases, const int offset)
{
    Lanes lanes;
    for (int l = 0; l < L; ++l) {
        lanes[l] = cases[(size_t(l) + size_t(offset)) % cases.size()];
    }
    return lanes;
}

// ---------------------------------------------------------------------------
// Comparison

/// @brief The bound one lane must meet, as a Catch2 `Approx`.
///
/// @warning `epsilon(0)` is not optional. `Approx` **ors** its margin test
/// with an epsilon test, and its default epsilon is float-grade
/// (`FLT_EPSILON * 100`, about 1.19e-5 relative), so a bare
/// `Approx(x).margin(m)` would also accept anything within 1.19e-5 of `x` --
/// silently defeating a 1e-14 double comparison. Zeroing it leaves the
/// absolute `margin` as the only bound.
///
/// `Approx` compares without subtracting, so `±inf` matches only itself and
/// only with the same sign, and a NaN never compares equal to anything. That
/// is the contract the barrier tests rely on: a barrier at penetration must
/// resolve to `+inf` or 0, never NaN.
inline Catch::Approx
approx(const double expected, const double tol, const double scale)
{
    // A non-finite scale must not reach the margin. `tol * inf` is `inf`, and
    // an infinite margin compares equal to *everything* -- including a finite
    // value where `+inf` was expected, which is precisely the case the barrier
    // tests rely on catching.
    const double margin = std::isfinite(scale) ? tol * scale : tol;
    return Catch::Approx(expected).epsilon(0).margin(margin);
}

/// @brief `approx` relative to the magnitude of `expected` itself.
inline Catch::Approx approx(const double expected, const double tol = VALUE_TOL)
{
    return approx(expected, tol, std::max(1.0, std::abs(expected)));
}

/// @brief Compare a batch-valued Eigen object against the scalar result, lane
/// by lane, where `scalar_of(l)` produces lane `l`'s scalar result.
template <typename BatchMatrix, typename ScalarOf>
void check_lanes(
    const std::string& name,
    const BatchMatrix& batched,
    ScalarOf&& scalar_of,
    const double tol = DERIVATIVE_TOL)
{
    for (int l = 0; l < L; ++l) {
        const auto expected = scalar_of(l);
        REQUIRE(batched.rows() == expected.rows());
        REQUIRE(batched.cols() == expected.cols());
        const double scale = std::max(1.0, expected.array().abs().maxCoeff());
        for (Eigen::Index i = 0; i < expected.size(); ++i) {
            CAPTURE(name, l, i, scale, tol);
            CHECK(batched(i).get(l) == approx(expected(i), tol, scale));
        }
    }
}

/// @brief Compare a batch of scalars against the scalar result, lane by lane.
template <typename ScalarOf>
void check_scalar_lanes(
    const std::string& name,
    const Batch& batched,
    ScalarOf&& scalar_of,
    const double tol = VALUE_TOL)
{
    const Lanes actual = Lanes::from(batched);
    for (int l = 0; l < L; ++l) {
        CAPTURE(name, l, tol);
        CHECK(actual[l] == approx(scalar_of(l), tol));
    }
}

/// @brief Sweep `cases` across every lane offset, comparing a one-argument
/// batch call against its scalar counterpart on each lane.
template <typename Container, typename ScalarFn, typename BatchFn>
void check_swept_lanes(
    const std::string& name,
    const Container& cases,
    ScalarFn&& scalar_fn,
    BatchFn&& batch_fn,
    const double tol = VALUE_TOL)
{
    for (int offset = 0; offset < int(cases.size()); ++offset) {
        const Lanes xs = lane_cases(cases, offset);
        const Lanes actual = Lanes::from(batch_fn(xs.load()));
        for (int l = 0; l < L; ++l) {
            CAPTURE(name, offset, l, xs[l], tol);
            CHECK(actual[l] == approx(scalar_fn(xs[l]), tol));
        }
    }
}

} // namespace ipc::tests

#endif
