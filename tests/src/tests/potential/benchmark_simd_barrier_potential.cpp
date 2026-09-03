// Per-collision barrier potential (value, gradient, Hessian) evaluated with
// four scalar types -- double, float, SimdBatch<double>, SimdBatch<float> --
// on the same collision set, so the four are compared on identical work. Runs
// on every scene in `assembly_scene_specs()` that is available, except the two
// smallest.
//
// The chain each variant evaluates is the one `NormalPotential` uses:
//
//   d = dist²(x),  b = κ·w·f(d),  ∇b = κ·w·f'(d)·∇d,
//   ∇²b = κ·w·(f"(d)·∇d∇dᵀ + f'(d)·∇²d)
//
// with the distance type resolved scalar-side (from the collision set) and
// passed explicitly, which a batch requires. Collisions are grouped by kind
// and distance type, then packed into an array-of-structures-of-arrays layout
// with one collision per lane, so a batch load is one contiguous read.
//
// The edge-edge mollifier is *not* applied by any variant: its
// implementation branches on a scalar `if` and is only instantiated for
// `float`/`double`, so a batch cannot evaluate it. The scene report below
// says how many edge-edge collisions are actually mollified (m < 1); for the
// rest the mollified and unmollified derivatives are identical, which is what
// the check against the library path relies on.
//
// The float variants come in two flavours. "float" converts the scene's
// coordinates as they are. "float (rescaled)" first re-centers each stencil
// on its centroid and divides by d̂, in double, so a float sees O(1) numbers:
// the distance functions are translation invariant and the barrier is
// homogeneous in (d, d̂), so the result is recovered exactly by scaling κ. The
// two flavours run the very same instructions; only the packing differs.
//
// Run with (single-threaded is the primary measurement):
//   ./ipc_toolkit_tests "[simd_barrier_potential]"
//
// Environment:
//   IPC_TOOLKIT_BENCH_SAMPLES  number of timed runs per cell (default 5)
//   IPC_TOOLKIT_BENCH_OUTPUT   write the results as JSON to this path (one
//                              report per scene under "scenes")

#include "assembly_scene.hpp"

#include <catch2/catch_test_macros.hpp>

#include <ipc/barrier/barrier.hpp>
#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/edge_edge_mollifier.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/distance/point_triangle.hpp>
#include <ipc/utils/simd.hpp>

#include <tbb/blocked_range.h>
#include <tbb/parallel_reduce.h>
#include <tbb/task_arena.h>

#include <spdlog/fmt/fmt.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

using namespace ipc;

namespace {

// -- Scalar-type abstraction --------------------------------------------------

/// @brief What the evaluation loop needs to know about a scalar type `T`.
template <typename T> struct ScalarTraits {
    using Real = T;
    static constexpr int LANES = 1;
    static T load(const Real* p) { return *p; }
    static void store(Real* p, const T& v) { *p = v; }
    static double reduce_add(const T& v) { return double(v); }
};

#ifdef IPC_TOOLKIT_WITH_SIMD
template <typename R, typename A> struct ScalarTraits<xsimd::batch<R, A>> {
    using Batch = xsimd::batch<R, A>;
    using Real = R;
    static constexpr int LANES = int(Batch::size);
    static Batch load(const Real* p) { return Batch::load_unaligned(p); }
    static void store(Real* p, const Batch& v) { v.store_unaligned(p); }
    static double reduce_add(const Batch& v)
    {
        return double(xsimd::reduce_add(v));
    }
};
#endif

template <typename T> std::string variant_name(const bool rescaled)
{
    using Tr = ScalarTraits<T>;
    const bool is_float = std::is_same_v<typename Tr::Real, float>;
    std::string name;
    if constexpr (Tr::LANES == 1) {
        name = is_float ? "float" : "double";
    } else {
        name = fmt::format(
            "simd<{}> x{}", is_float ? "float" : "double", Tr::LANES);
    }
    return rescaled ? name + " (rescaled)" : name;
}

// -- Collision kinds ----------------------------------------------------------

enum class Kind : uint8_t { VV, EV, EE, FV };

constexpr int num_vertices(const Kind k)
{
    switch (k) {
    case Kind::VV:
        return 2;
    case Kind::EV:
        return 3;
    default:
        return 4;
    }
}

const char* kind_name(const Kind k)
{
    switch (k) {
    case Kind::VV:
        return "vertex-vertex";
    case Kind::EV:
        return "edge-vertex";
    case Kind::EE:
        return "edge-edge";
    default:
        return "face-vertex";
    }
}

/// @brief The distance query for a kind, for any scalar `T`.
template <typename T, Kind K> struct Stencil;

template <typename T> struct Stencil<T, Kind::VV> {
    static constexpr int NV = 2;
    using V = Eigen::Vector3<T>;
    static T distance(const V* x, int /*dtype*/)
    {
        return point_point_distance(x[0], x[1]);
    }
    static Eigen::Vector<T, 6> gradient(const V* x, int /*dtype*/)
    {
        return point_point_distance_gradient(x[0], x[1]);
    }
    static Eigen::Matrix<T, 6, 6> hessian(const V* x, int /*dtype*/)
    {
        return point_point_distance_hessian(x[0], x[1]);
    }
};

template <typename T> struct Stencil<T, Kind::EV> {
    static constexpr int NV = 3;
    using V = Eigen::Vector3<T>;
    static T distance(const V* x, int dtype)
    {
        return point_edge_distance(
            x[0], x[1], x[2], PointEdgeDistanceType(dtype));
    }
    static Eigen::Vector<T, 9> gradient(const V* x, int dtype)
    {
        return point_edge_distance_gradient(
            x[0], x[1], x[2], PointEdgeDistanceType(dtype));
    }
    static Eigen::Matrix<T, 9, 9> hessian(const V* x, int dtype)
    {
        return point_edge_distance_hessian(
            x[0], x[1], x[2], PointEdgeDistanceType(dtype));
    }
};

template <typename T> struct Stencil<T, Kind::EE> {
    static constexpr int NV = 4;
    using V = Eigen::Vector3<T>;
    static T distance(const V* x, int dtype)
    {
        return edge_edge_distance(
            x[0], x[1], x[2], x[3], EdgeEdgeDistanceType(dtype));
    }
    static Eigen::Vector<T, 12> gradient(const V* x, int dtype)
    {
        return edge_edge_distance_gradient(
            x[0], x[1], x[2], x[3], EdgeEdgeDistanceType(dtype));
    }
    static Eigen::Matrix<T, 12, 12> hessian(const V* x, int dtype)
    {
        return edge_edge_distance_hessian(
            x[0], x[1], x[2], x[3], EdgeEdgeDistanceType(dtype));
    }
};

template <typename T> struct Stencil<T, Kind::FV> {
    static constexpr int NV = 4;
    using V = Eigen::Vector3<T>;
    static T distance(const V* x, int dtype)
    {
        return point_triangle_distance(
            x[0], x[1], x[2], x[3], PointTriangleDistanceType(dtype));
    }
    static Eigen::Vector<T, 12> gradient(const V* x, int dtype)
    {
        return point_triangle_distance_gradient(
            x[0], x[1], x[2], x[3], PointTriangleDistanceType(dtype));
    }
    static Eigen::Matrix<T, 12, 12> hessian(const V* x, int dtype)
    {
        return point_triangle_distance_hessian(
            x[0], x[1], x[2], x[3], PointTriangleDistanceType(dtype));
    }
};

// -- Grouping and packing -----------------------------------------------------

/// @brief The collisions of one kind and distance type (lane-layout agnostic).
struct GroupSpec {
    Kind kind;
    int dtype;
    std::vector<size_t> ids; ///< Indices into the NormalCollisions.

    std::string name() const
    {
        return fmt::format("{} (dtype {})", kind_name(kind), dtype);
    }
};

/// @brief One group packed for a lane count `L`, in real type `R`.
///
/// Layout: block b (L consecutive collisions) occupies
/// `x[b*3*NV*L + comp*L + lane]`, i.e. within a block each of the 3·NV
/// coordinates is stored contiguously across lanes. For L = 1 this is the
/// plain per-collision DOF vector `CollisionStencil::dof` returns.
template <typename R> struct PackedGroup {
    Kind kind = Kind::VV;
    int dtype = 0;
    size_t n = 0;       ///< Real collisions.
    size_t nblocks = 0; ///< ceil(n / L).
    int lanes = 1;
    std::vector<R> x; ///< nblocks * 3*NV * L
    std::vector<R> w; ///< nblocks * L; zero on padding lanes.
    std::vector<R> grad_out;
    std::vector<R> hess_out;

    int nv() const { return num_vertices(kind); }
    int ndof() const { return 3 * nv(); }

    R grad_at(size_t i, int comp) const
    {
        return grad_out
            [(i / lanes) * ndof() * lanes + comp * lanes + (i % lanes)];
    }
    R hess_at(size_t i, int entry) const
    {
        return hess_out
            [(i / lanes) * ndof() * ndof() * lanes + entry * lanes
             + (i % lanes)];
    }
};

std::vector<GroupSpec> group_collisions(const NormalCollisions& collisions)
{
    // A map from (kind, dtype) to its group, kept in first-seen order.
    std::vector<GroupSpec> groups;
    const auto group_for = [&](const Kind kind, const int dtype) -> GroupSpec& {
        for (GroupSpec& g : groups) {
            if (g.kind == kind && g.dtype == dtype) {
                return g;
            }
        }
        groups.push_back(GroupSpec { kind, dtype, {} });
        return groups.back();
    };

    for (size_t i = 0; i < collisions.size(); i++) {
        if (collisions.is_vertex_vertex(i)) {
            group_for(Kind::VV, int(PointPointDistanceType::P_P))
                .ids.push_back(i);
        } else if (collisions.is_edge_vertex(i)) {
            const auto& ev =
                dynamic_cast<const EdgeVertexNormalCollision&>(collisions[i]);
            group_for(Kind::EV, int(ev.known_dtype())).ids.push_back(i);
        } else if (collisions.is_edge_edge(i)) {
            const auto& ee =
                dynamic_cast<const EdgeEdgeNormalCollision&>(collisions[i]);
            group_for(Kind::EE, int(ee.known_dtype())).ids.push_back(i);
        } else if (collisions.is_face_vertex(i)) {
            const auto& fv =
                dynamic_cast<const FaceVertexNormalCollision&>(collisions[i]);
            group_for(Kind::FV, int(fv.known_dtype())).ids.push_back(i);
        } else {
            throw std::runtime_error("unsupported collision kind");
        }
    }

    // Largest groups first, so a truncated run still covers the bulk.
    std::stable_sort(
        groups.begin(), groups.end(),
        [](const GroupSpec& a, const GroupSpec& b) {
            return a.ids.size() > b.ids.size();
        });
    return groups;
}

/// @param scale Coordinates are stored as (x - centroid) / scale; 1 means as
///              they are (no re-centering).
template <typename R>
PackedGroup<R> pack_group(
    const GroupSpec& spec,
    const int lanes,
    const ipc::tests::AssemblyScene& scene,
    const double scale)
{
    const NormalCollisions& collisions = scene.collisions();
    const Eigen::MatrixXi& edges = scene.mesh().edges();
    const Eigen::MatrixXi& faces = scene.mesh().faces();
    const Eigen::MatrixXd& X = scene.vertices();

    PackedGroup<R> g;
    g.kind = spec.kind;
    g.dtype = spec.dtype;
    g.lanes = lanes;
    g.n = spec.ids.size();
    g.nblocks = (g.n + lanes - 1) / lanes;
    const int ndof = g.ndof();

    g.x.assign(g.nblocks * ndof * lanes, R(0));
    g.w.assign(g.nblocks * lanes, R(0));

    for (size_t i = 0; i < g.nblocks * size_t(lanes); i++) {
        // Padding lanes replay the last collision with zero weight so every
        // lane holds a finite, in-range configuration.
        const size_t src = std::min(i, g.n - 1);
        const NormalCollision& c = collisions[spec.ids[src]];
        VectorMax12d dof = c.dof(X, edges, faces);
        assert(dof.size() == ndof);

        if (scale != 1.0) {
            const Eigen::Vector3d centroid =
                dof.reshaped(3, dof.size() / 3).rowwise().mean();
            for (int v = 0; v < ndof / 3; v++) {
                dof.segment<3>(3 * v) =
                    (dof.segment<3>(3 * v) - centroid) / scale;
            }
        }

        R* block = g.x.data() + (i / lanes) * ndof * lanes + (i % lanes);
        for (int k = 0; k < ndof; k++) {
            block[k * lanes] = R(dof[k]);
        }
        g.w[i] = i < g.n ? R(c.weight) : R(0);
    }
    return g;
}

// -- The evaluation loop ------------------------------------------------------

struct Params {
    double dhat_sqr; ///< The potential is a function of squared distance.
    double kappa;
};

enum class Quantity : uint8_t { VALUE, GRADIENT, HESSIAN, HESSIAN_SUM };

const char* quantity_name(const Quantity q)
{
    switch (q) {
    case Quantity::VALUE:
        return "value";
    case Quantity::GRADIENT:
        return "gradient";
    case Quantity::HESSIAN:
        return "hessian";
    default:
        return "hessian_sum";
    }
}

constexpr Quantity QUANTITIES[] = { Quantity::VALUE, Quantity::GRADIENT,
                                    Quantity::HESSIAN, Quantity::HESSIAN_SUM };

template <typename T, Kind K> struct Evaluator {
    using Tr = ScalarTraits<T>;
    using R = typename Tr::Real;
    using S = Stencil<T, K>;
    static constexpr int L = Tr::LANES;
    static constexpr int NV = S::NV;
    static constexpr int N = 3 * NV;
    /// @brief Blocks summed in T before the partial sum is moved to double.
    static constexpr size_t REDUCE_CHUNK = 256;

    static void
    load(const PackedGroup<R>& g, const size_t b, Eigen::Vector3<T>* x, T& w)
    {
        const R* in = g.x.data() + b * (N * L);
        for (int v = 0; v < NV; v++) {
            for (int k = 0; k < 3; k++) {
                x[v][k] = Tr::load(in + (3 * v + k) * L);
            }
        }
        w = Tr::load(g.w.data() + b * L);
    }

    static double
    value(const PackedGroup<R>& g, const Params& p, size_t b0, size_t b1)
    {
        const T dhat(R(p.dhat_sqr));
        const T kappa(R(p.kappa));
        Eigen::Vector3<T> x[NV];
        T w;
        double total = 0;
        // Sum in T over a chunk, then in double across chunks, so a float's
        // accumulation error is bounded by the chunk length and the accuracy
        // column reflects the per-collision precision.
        for (size_t c0 = b0; c0 < b1; c0 += REDUCE_CHUNK) {
            T acc(R(0));
            const size_t c1 = std::min(b1, c0 + REDUCE_CHUNK);
            for (size_t b = c0; b < c1; b++) {
                load(g, b, x, w);
                const T d = S::distance(x, g.dtype);
                acc += (kappa * w) * barrier(d, dhat);
            }
            total += Tr::reduce_add(acc);
        }
        return total;
    }

    static void
    gradient(PackedGroup<R>& g, const Params& p, size_t b0, size_t b1)
    {
        const T dhat(R(p.dhat_sqr));
        const T kappa(R(p.kappa));
        Eigen::Vector3<T> x[NV];
        T w;
        for (size_t b = b0; b < b1; b++) {
            load(g, b, x, w);
            const T d = S::distance(x, g.dtype);
            const Eigen::Vector<T, N> grad_d = S::gradient(x, g.dtype);
            const T grad_f = barrier_first_derivative(d, dhat);
            const Eigen::Vector<T, N> grad = (kappa * w * grad_f) * grad_d;

            R* out = g.grad_out.data() + b * (N * L);
            for (int k = 0; k < N; k++) {
                Tr::store(out + k * L, grad[k]);
            }
        }
    }

    static Eigen::Matrix<T, N, N> local_hessian(
        const Eigen::Vector3<T>* x,
        const T& w,
        const int dtype,
        const T& dhat,
        const T& kappa)
    {
        const T d = S::distance(x, dtype);
        const Eigen::Vector<T, N> grad_d = S::gradient(x, dtype);
        const Eigen::Matrix<T, N, N> hess_d = S::hessian(x, dtype);
        const T grad_f = barrier_first_derivative(d, dhat);
        const T hess_f = barrier_second_derivative(d, dhat);
        const T kw = kappa * w;
        return (kw * hess_f) * (grad_d * grad_d.transpose())
            + (kw * grad_f) * hess_d;
    }

    static void
    hessian(PackedGroup<R>& g, const Params& p, size_t b0, size_t b1)
    {
        const T dhat(R(p.dhat_sqr));
        const T kappa(R(p.kappa));
        Eigen::Vector3<T> x[NV];
        T w;
        for (size_t b = b0; b < b1; b++) {
            load(g, b, x, w);
            const Eigen::Matrix<T, N, N> hess =
                local_hessian(x, w, g.dtype, dhat, kappa);
            R* out = g.hess_out.data() + b * (N * N * L);
            for (int k = 0; k < N * N; k++) {
                Tr::store(out + k * L, hess.data()[k]);
            }
        }
    }

    /// @brief The Hessian's compute cost alone: every entry is summed into
    /// one accumulator instead of being written out, so the 144 stores per
    /// collision -- which are memory-bound once threads share the bus -- are
    /// not part of the measurement.
    static double
    hessian_sum(const PackedGroup<R>& g, const Params& p, size_t b0, size_t b1)
    {
        const T dhat(R(p.dhat_sqr));
        const T kappa(R(p.kappa));
        Eigen::Vector3<T> x[NV];
        T w;
        double total = 0;
        for (size_t c0 = b0; c0 < b1; c0 += REDUCE_CHUNK) {
            T acc(R(0));
            const size_t c1 = std::min(b1, c0 + REDUCE_CHUNK);
            for (size_t b = c0; b < c1; b++) {
                load(g, b, x, w);
                acc += local_hessian(x, w, g.dtype, dhat, kappa).sum();
            }
            total += Tr::reduce_add(acc);
        }
        return total;
    }
};

/// @brief Dispatch on the kind at runtime; everything inside is static.
template <typename T>
double run_group(
    PackedGroup<typename ScalarTraits<T>::Real>& g,
    const Params& p,
    const Quantity q,
    const size_t b0,
    const size_t b1)
{
    const auto run = [&](auto kind_tag) -> double {
        using E = Evaluator<T, decltype(kind_tag)::value>;
        switch (q) {
        case Quantity::VALUE:
            return E::value(g, p, b0, b1);
        case Quantity::GRADIENT:
            E::gradient(g, p, b0, b1);
            return 0;
        case Quantity::HESSIAN:
            E::hessian(g, p, b0, b1);
            return 0;
        default:
            return E::hessian_sum(g, p, b0, b1);
        }
    };
    switch (g.kind) {
    case Kind::VV:
        return run(std::integral_constant<Kind, Kind::VV>());
    case Kind::EV:
        return run(std::integral_constant<Kind, Kind::EV>());
    case Kind::EE:
        return run(std::integral_constant<Kind, Kind::EE>());
    default:
        return run(std::integral_constant<Kind, Kind::FV>());
    }
}

/// @brief All groups packed for one scalar type, plus its outputs.
template <typename T> struct Variant {
    using Tr = ScalarTraits<T>;
    using R = typename Tr::Real;

    /// @brief Coordinates are stored as (x - centroid) / scale.
    double scale = 1.0;
    std::vector<PackedGroup<R>> groups;

    void pack(
        const std::vector<GroupSpec>& specs,
        const ipc::tests::AssemblyScene& scene)
    {
        groups.clear();
        for (const GroupSpec& spec : specs) {
            groups.push_back(pack_group<R>(spec, Tr::LANES, scene, scale));
        }
    }

    void allocate_outputs(const Quantity q)
    {
        for (PackedGroup<R>& g : groups) {
            const size_t per_block = size_t(g.ndof()) * g.lanes;
            if (q == Quantity::GRADIENT) {
                g.grad_out.assign(g.nblocks * per_block, R(0));
            } else if (q == Quantity::HESSIAN) {
                g.hess_out.assign(g.nblocks * per_block * g.ndof(), R(0));
            }
        }
    }

    /// @brief The parameters in the packed coordinates.
    ///
    /// With x' = (x - c)/s: d'² = d²/s², and since the barrier is homogeneous,
    /// b(d'², d̂²/s²) = b(d², d̂²)/s⁴. Its first derivative wrt d² picks up a
    /// further s², ∇ₓ' a further s, so the value, gradient, and Hessian are
    /// s⁴, s³, and s² too small; folding those into κ recovers them exactly.
    Params effective(const Params& p, const Quantity q) const
    {
        if (scale == 1.0) {
            return p;
        }
        const double s = scale;
        const double kappa_scale = q == Quantity::VALUE ? s * s * s * s
            : q == Quantity::GRADIENT                   ? s * s * s
                                                        : s * s;
        return Params { p.dhat_sqr / (s * s), p.kappa * kappa_scale };
    }

    double run(const Params& params, const Quantity q, const bool parallel)
    {
        const Params p = effective(params, q);
        double total = 0;
        for (PackedGroup<R>& g : groups) {
            if (!parallel) {
                total += run_group<T>(g, p, q, 0, g.nblocks);
                continue;
            }
            // Roughly 1k collisions per task, so the scheduling overhead is
            // negligible next to the work.
            const size_t grain = std::max<size_t>(1, 1024 / Tr::LANES);
            total += tbb::parallel_reduce(
                tbb::blocked_range<size_t>(0, g.nblocks, grain), 0.0,
                [&](const tbb::blocked_range<size_t>& r, double partial) {
                    return partial + run_group<T>(g, p, q, r.begin(), r.end());
                },
                std::plus<double>());
        }
        return total;
    }
};

// -- Timing -------------------------------------------------------------------

struct Timing {
    double median_s = 0;
    double min_s = 0;
};

template <typename F> Timing time_runs(F&& f, const int num_samples)
{
    std::vector<double> samples;
    samples.reserve(num_samples);
    for (int i = 0; i < num_samples; i++) {
        const auto start = std::chrono::steady_clock::now();
        f();
        const auto end = std::chrono::steady_clock::now();
        samples.push_back(std::chrono::duration<double>(end - start).count());
    }
    std::sort(samples.begin(), samples.end());
    return Timing { samples[samples.size() / 2], samples.front() };
}

/// @brief |got - want| / |want|, where two infinities of the same sign agree:
/// a scene with a zero-distance collision has an infinite barrier value.
double relative_error(const double got, const double want)
{
    if (!std::isfinite(want)) {
        return got == want ? 0.0 : std::numeric_limits<double>::infinity();
    }
    return std::abs(got - want) / std::abs(want);
}

int env_int(const char* name, const int fallback)
{
    const char* s = std::getenv(name);
    return s != nullptr ? std::atoi(s) : fallback;
}

// -- Accuracy -----------------------------------------------------------------

struct Accuracy {
    double value_rel = 0;
    double grad_max_rel = 0;
    double grad_median_rel = 0;
    double hess_max_rel = 0;
    double hess_median_rel = 0;
    size_t non_finite = 0;
    /// @brief Collisions with at least one non-finite entry, per group.
    std::vector<std::pair<std::string, size_t>> non_finite_by_group;
};

double median_of(std::vector<double>& v)
{
    if (v.empty()) {
        return 0;
    }
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
}

/// @brief Per-collision relative error in the ∞-norm, scaled by the
/// reference's ∞-norm (derivatives of near-parallel edges are ill-conditioned,
/// so an entry-wise relative error would be dominated by cancellations).
/// @return Number of collisions with a non-finite entry.
template <typename R, typename RRef>
size_t compare_outputs(
    const PackedGroup<R>& got,
    const PackedGroup<RRef>& ref,
    const Quantity q,
    std::vector<double>& rel_errors,
    size_t& non_finite)
{
    const int n_entries =
        q == Quantity::GRADIENT ? got.ndof() : got.ndof() * got.ndof();
    size_t bad_collisions = 0;
    for (size_t i = 0; i < got.n; i++) {
        double max_diff = 0, max_ref = 0;
        bool any_non_finite = false;
        for (int k = 0; k < n_entries; k++) {
            const double a = q == Quantity::GRADIENT
                ? double(got.grad_at(i, k))
                : double(got.hess_at(i, k));
            const double b = q == Quantity::GRADIENT
                ? double(ref.grad_at(i, k))
                : double(ref.hess_at(i, k));
            if (!std::isfinite(a)) {
                non_finite++;
                any_non_finite = true;
            }
            max_diff = std::max(max_diff, std::abs(a - b));
            max_ref = std::max(max_ref, std::abs(b));
        }
        bad_collisions += any_non_finite;
        rel_errors.push_back(
            max_diff / std::max(max_ref, std::numeric_limits<double>::min()));
    }
    return bad_collisions;
}

// -- The library path, for reference ------------------------------------------

/// @brief The library's own per-collision evaluation (virtual dispatch,
/// dynamic-size `VectorMax12d`, and the mollifier for edge-edge), gathering
/// the DOF from the vertex matrix as `Potential::gradient` does.
double library_path(
    const ipc::tests::AssemblyScene& scene,
    const Quantity q,
    const bool parallel,
    std::vector<double>& out)
{
    const NormalCollisions& collisions = scene.collisions();
    const BarrierPotential& potential = scene.potential();
    const Eigen::MatrixXi& edges = scene.mesh().edges();
    const Eigen::MatrixXi& faces = scene.mesh().faces();
    const Eigen::MatrixXd& X = scene.vertices();

    const size_t stride = q == Quantity::GRADIENT ? 12 : 144;
    const bool stores = q == Quantity::GRADIENT || q == Quantity::HESSIAN;
    if (stores && out.size() != collisions.size() * stride) {
        out.assign(collisions.size() * stride, 0.0);
    }

    const auto body = [&](const size_t begin, const size_t end) {
        double partial = 0;
        for (size_t i = begin; i < end; i++) {
            const NormalCollision& c = collisions[i];
            const VectorMax12d dof = c.dof(X, edges, faces);
            if (q == Quantity::VALUE) {
                partial += potential(c, dof);
            } else if (q == Quantity::GRADIENT) {
                const VectorMax12d grad = potential.gradient(c, dof);
                std::copy(
                    grad.data(), grad.data() + grad.size(),
                    out.data() + i * stride);
            } else if (q == Quantity::HESSIAN) {
                const MatrixMax12d hess =
                    potential.hessian(c, dof, PSDProjectionMethod::NONE);
                std::copy(
                    hess.data(), hess.data() + hess.size(),
                    out.data() + i * stride);
            } else {
                partial +=
                    potential.hessian(c, dof, PSDProjectionMethod::NONE).sum();
            }
        }
        return partial;
    };

    if (!parallel) {
        return body(0, collisions.size());
    }
    return tbb::parallel_reduce(
        tbb::blocked_range<size_t>(0, collisions.size(), 1024), 0.0,
        [&](const tbb::blocked_range<size_t>& r, double partial) {
            return partial + body(r.begin(), r.end());
        },
        std::plus<double>());
}

// -- Results ------------------------------------------------------------------

/// @brief A double as a JSON value; JSON has no inf/nan, so those are strings.
std::string json_number(const double x)
{
    return std::isfinite(x) ? fmt::format("{:.6g}", x)
                            : fmt::format("\"{}\"", x);
}

struct Cell {
    std::string variant;
    Quantity quantity;
    bool parallel;
    Timing timing;
};

struct Report {
    std::string scene;
    size_t num_collisions = 0;
    std::vector<std::pair<std::string, size_t>> composition;
    size_t num_mollified_ee = 0;
    int num_samples = 0;
    int num_threads = 1;
    std::vector<Cell> cells;
    std::vector<std::pair<std::string, double>> pack_seconds;
    std::vector<std::pair<std::string, Accuracy>> accuracy;
    std::string reference_check;

    std::string to_json() const
    {
        fmt::memory_buffer buf;
        auto f = std::back_inserter(buf);
        fmt::format_to(f, "{{\n");
        fmt::format_to(f, "  \"scene\": \"{}\",\n", scene);
        fmt::format_to(f, "  \"num_collisions\": {},\n", num_collisions);
        fmt::format_to(f, "  \"num_mollified_ee\": {},\n", num_mollified_ee);
        fmt::format_to(f, "  \"num_samples\": {},\n", num_samples);
        fmt::format_to(f, "  \"num_threads\": {},\n", num_threads);
        fmt::format_to(f, "  \"composition\": {{");
        for (size_t i = 0; i < composition.size(); i++) {
            fmt::format_to(
                f, "{}\"{}\": {}", i ? ", " : "", composition[i].first,
                composition[i].second);
        }
        fmt::format_to(f, "}},\n  \"pack_seconds\": {{");
        for (size_t i = 0; i < pack_seconds.size(); i++) {
            fmt::format_to(
                f, "{}\"{}\": {:.9g}", i ? ", " : "", pack_seconds[i].first,
                pack_seconds[i].second);
        }
        fmt::format_to(f, "}},\n  \"accuracy\": {{\n");
        for (size_t i = 0; i < accuracy.size(); i++) {
            const Accuracy& a = accuracy[i].second;
            fmt::format_to(
                f,
                "    \"{}\": {{\"value_rel\": {}, \"grad_max_rel\": {}, "
                "\"grad_median_rel\": {}, \"hess_max_rel\": {}, "
                "\"hess_median_rel\": {}, \"non_finite\": {}, "
                "\"non_finite_collisions_by_group\": {{",
                accuracy[i].first, json_number(a.value_rel),
                json_number(a.grad_max_rel), json_number(a.grad_median_rel),
                json_number(a.hess_max_rel), json_number(a.hess_median_rel),
                a.non_finite);
            for (size_t j = 0; j < a.non_finite_by_group.size(); j++) {
                fmt::format_to(
                    f, "{}\"{}\": {}", j ? ", " : "",
                    a.non_finite_by_group[j].first,
                    a.non_finite_by_group[j].second);
            }
            fmt::format_to(f, "}}}}{}\n", i + 1 < accuracy.size() ? "," : "");
        }
        fmt::format_to(f, "  }},\n  \"timings\": [\n");
        for (size_t i = 0; i < cells.size(); i++) {
            const Cell& c = cells[i];
            fmt::format_to(
                f,
                "    {{\"variant\": \"{}\", \"quantity\": \"{}\", "
                "\"parallel\": {}, \"median_s\": {:.9g}, \"min_s\": {:.9g}}}{}"
                "\n",
                c.variant, quantity_name(c.quantity),
                c.parallel ? "true" : "false", c.timing.median_s,
                c.timing.min_s, i + 1 < cells.size() ? "," : "");
        }
        fmt::format_to(f, "  ]\n}}");
        return fmt::to_string(buf);
    }
};

/// @brief Time and check one scalar type; append to the report.
template <typename T>
void bench_variant(
    const ipc::tests::AssemblyScene& scene,
    const std::vector<GroupSpec>& specs,
    const Params& params,
    Variant<double>& reference,
    Report& report,
    Variant<T>& variant)
{
    const std::string name = variant_name<T>(variant.scale != 1.0);

    const Timing pack_time = time_runs([&] { variant.pack(specs, scene); }, 1);
    report.pack_seconds.emplace_back(name, pack_time.median_s);

    Accuracy acc;
    std::vector<size_t> bad_by_group(specs.size(), 0);
    for (const bool parallel : { false, true }) {
        for (const Quantity q : QUANTITIES) {
            variant.allocate_outputs(q);
            // Warm-up touches the freshly allocated outputs so first-touch
            // page faults are not in the timed runs.
            double value = variant.run(params, q, parallel);
            const Timing t = time_runs(
                [&] { value = variant.run(params, q, parallel); },
                report.num_samples);
            report.cells.push_back(Cell { name, q, parallel, t });

            if (parallel) {
                continue;
            }
            if (q == Quantity::VALUE) {
                const double ref_value = reference.run(params, q, false);
                acc.value_rel = relative_error(value, ref_value);
            } else if (q != Quantity::HESSIAN_SUM) {
                std::vector<double> rel;
                for (size_t i = 0; i < variant.groups.size(); i++) {
                    bad_by_group[i] += compare_outputs(
                        variant.groups[i], reference.groups[i], q, rel,
                        acc.non_finite);
                }
                const double max_rel =
                    rel.empty() ? 0 : *std::max_element(rel.begin(), rel.end());
                const double median_rel = median_of(rel);
                if (q == Quantity::GRADIENT) {
                    acc.grad_max_rel = max_rel;
                    acc.grad_median_rel = median_rel;
                } else {
                    acc.hess_max_rel = max_rel;
                    acc.hess_median_rel = median_rel;
                }
            }
        }
    }
    for (size_t i = 0; i < specs.size(); i++) {
        if (bad_by_group[i] > 0) {
            acc.non_finite_by_group.emplace_back(
                specs[i].name(), bad_by_group[i]);
        }
    }
    report.accuracy.emplace_back(name, acc);
}

/// @brief Check the double path against the library on every collision the
/// mollifier leaves untouched (m = 1), where the two must agree to rounding.
std::string check_against_library(
    const ipc::tests::AssemblyScene& scene,
    const std::vector<GroupSpec>& specs,
    const Params& params,
    Variant<double>& variant,
    size_t& num_mollified_ee)
{
    const NormalCollisions& collisions = scene.collisions();
    const BarrierPotential& potential = scene.potential();
    const Eigen::MatrixXi& edges = scene.mesh().edges();
    const Eigen::MatrixXi& faces = scene.mesh().faces();
    const Eigen::MatrixXd& X = scene.vertices();

    variant.allocate_outputs(Quantity::GRADIENT);
    variant.run(params, Quantity::GRADIENT, true);
    variant.allocate_outputs(Quantity::HESSIAN);
    variant.run(params, Quantity::HESSIAN, true);
    const double value = variant.run(params, Quantity::VALUE, true);

    double ref_value = 0, max_grad_rel = 0, max_hess_rel = 0;
    size_t checked = 0;
    num_mollified_ee = 0;
    for (size_t gi = 0; gi < specs.size(); gi++) {
        const GroupSpec& spec = specs[gi];
        const PackedGroup<double>& g = variant.groups[gi];
        for (size_t i = 0; i < spec.ids.size(); i++) {
            const NormalCollision& c = collisions[spec.ids[i]];
            const VectorMax12d dof = c.dof(X, edges, faces);

            // The unmollified potential, from the library's own pieces.
            ref_value += params.kappa * c.weight
                * potential.barrier()(c.compute_distance(dof), params.dhat_sqr);

            if (c.is_mollified() && c.mollifier(dof) < 1.0) {
                num_mollified_ee++;
                continue;
            }

            const VectorMax12d grad = potential.gradient(c, dof);
            const MatrixMax12d hess =
                potential.hessian(c, dof, PSDProjectionMethod::NONE);
            double grad_diff = 0, hess_diff = 0;
            for (int k = 0; k < grad.size(); k++) {
                grad_diff =
                    std::max(grad_diff, std::abs(grad[k] - g.grad_at(i, k)));
            }
            for (int k = 0; k < hess.size(); k++) {
                hess_diff = std::max(
                    hess_diff, std::abs(hess.data()[k] - g.hess_at(i, k)));
            }
            const double grad_scale =
                std::max(grad.template lpNorm<Eigen::Infinity>(), 1e-300);
            const double hess_scale =
                std::max(hess.template lpNorm<Eigen::Infinity>(), 1e-300);
            max_grad_rel = std::max(max_grad_rel, grad_diff / grad_scale);
            max_hess_rel = std::max(max_hess_rel, hess_diff / hess_scale);
            checked++;
        }
    }

    return fmt::format(
        "double vs. library on {} unmollified collisions: value rel. err "
        "{:.2e}{}, gradient max rel. err {:.2e}, Hessian max rel. err {:.2e}",
        checked, relative_error(value, ref_value),
        std::isfinite(ref_value)
            ? ""
            : " (the value is infinite: a collision has zero distance)",
        max_grad_rel, max_hess_rel);
}

void print_report(const Report& r)
{
    fmt::print("\n=== Barrier potential per scalar type: {} ===\n", r.scene);
    fmt::print("{} collisions:", r.num_collisions);
    for (const auto& [name, count] : r.composition) {
        fmt::print(" {} {};", count, name);
    }
    fmt::print(
        "\n{} edge-edge collisions are mollified (m < 1) and evaluated here "
        "without the mollifier.\n",
        r.num_mollified_ee);
    fmt::print(
        "{}\nMedian of {} runs; parallel columns use {} threads. "
        "\"hess-sum\" evaluates the Hessian without storing it.\n\n",
        r.reference_check, r.num_samples, r.num_threads);

    const auto find = [&](const std::string& v, Quantity q, bool par) {
        for (const Cell& c : r.cells) {
            if (c.variant == v && c.quantity == q && c.parallel == par) {
                return c.timing.median_s;
            }
        }
        return std::numeric_limits<double>::quiet_NaN();
    };

    std::vector<std::string> variants;
    for (const Cell& c : r.cells) {
        if (std::find(variants.begin(), variants.end(), c.variant)
            == variants.end()) {
            variants.push_back(c.variant);
        }
    }

    const double per = 1e9 / double(r.num_collisions);
    for (const bool parallel : { false, true }) {
        fmt::print("--- {} ---\n", parallel ? "parallel" : "single-threaded");
        fmt::print("{:<24}", "variant");
        for (const char* h : { "value", "grad", "hess", "hess-sum" }) {
            fmt::print(" {:>9}", fmt::format("{}(ms)", h));
        }
        fmt::print(" |");
        for (const char* h : { "value", "grad", "hess", "hess-sum" }) {
            fmt::print(" {:>9}", fmt::format("{} ns", h));
        }
        fmt::print(" |");
        for (const char* h : { "value", "grad", "hess", "hess-sum" }) {
            fmt::print(" {:>9}", fmt::format("{} x", h));
        }
        fmt::print("\n");
        for (const std::string& v : variants) {
            fmt::print("{:<24}", v);
            for (const Quantity q : QUANTITIES) {
                fmt::print(" {:>9.2f}", find(v, q, parallel) * 1e3);
            }
            fmt::print(" |");
            for (const Quantity q : QUANTITIES) {
                fmt::print(" {:>9.1f}", find(v, q, parallel) * per);
            }
            fmt::print(" |");
            for (const Quantity q : QUANTITIES) {
                fmt::print(
                    " {:>8.2f}x",
                    find("double", q, parallel) / find(v, q, parallel));
            }
            fmt::print("\n");
        }
        fmt::print("\n");
    }

    fmt::print("--- packing (gather into lane layout, once per variant) ---\n");
    for (const auto& [name, s] : r.pack_seconds) {
        fmt::print("{:<24} {:>9.2f} ms\n", name, s * 1e3);
    }

    fmt::print("\n--- accuracy relative to double ---\n");
    fmt::print(
        "{:<24} {:>10} {:>12} {:>12} {:>12} {:>12} {:>10}\n", "variant",
        "value rel", "grad max", "grad median", "hess max", "hess median",
        "non-finite");
    for (const auto& [name, a] : r.accuracy) {
        fmt::print(
            "{:<24} {:>10.2e} {:>12.2e} {:>12.2e} {:>12.2e} {:>12.2e} {:>10}\n",
            name, a.value_rel, a.grad_max_rel, a.grad_median_rel,
            a.hess_max_rel, a.hess_median_rel, a.non_finite);
    }
    for (const auto& [name, a] : r.accuracy) {
        if (a.non_finite_by_group.empty()) {
            continue;
        }
        fmt::print("  {}: collisions with a non-finite entry:", name);
        for (const auto& [group, count] : a.non_finite_by_group) {
            fmt::print(" {} {};", count, group);
        }
        fmt::print("\n");
    }
    fmt::print("\n");
    std::fflush(stdout);
}

/// @brief Run every variant on one scene and return its report.
Report bench_scene(const ipc::tests::AssemblyScene& scene)
{
    Report report;
    report.scene = scene.label();
    report.num_collisions = scene.num_collisions();
    report.num_samples = env_int("IPC_TOOLKIT_BENCH_SAMPLES", 5);
    report.num_threads = tbb::this_task_arena::max_concurrency();

    const std::vector<GroupSpec> groups = group_collisions(scene.collisions());
    for (const GroupSpec& g : groups) {
        report.composition.emplace_back(g.name(), g.ids.size());
    }

    const double dhat = scene.potential().dhat();
    const Params params { dhat * dhat, scene.potential().stiffness() };

    // The double path first: it is the reference the others are checked
    // against, and is itself checked against the library.
    Variant<double> ref;
    ref.pack(groups, scene);
    report.reference_check = check_against_library(
        scene, groups, params, ref, report.num_mollified_ee);

    const auto bench = [&](auto& variant, const double scale) {
        variant.scale = scale;
        bench_variant(scene, groups, params, ref, report, variant);
    };
    {
        Variant<double> v;
        bench(v, 1.0);
    }
    {
        Variant<float> v;
        bench(v, 1.0);
        bench(v, dhat);
    }
#ifdef IPC_TOOLKIT_WITH_SIMD
    {
        Variant<SimdBatch<double>> v;
        bench(v, 1.0);
    }
    {
        Variant<SimdBatch<float>> v;
        bench(v, 1.0);
        bench(v, dhat);
    }
#endif

    // The library's own path, so the hand-rolled double loop above can be
    // placed relative to what `Potential::gradient`/`hessian` actually do.
    {
        std::vector<double> out;
        for (const bool parallel : { false, true }) {
            for (const Quantity q : QUANTITIES) {
                library_path(scene, q, parallel, out); // warm-up
                const Timing t = time_runs(
                    [&] { library_path(scene, q, parallel, out); },
                    report.num_samples);
                report.cells.push_back(
                    Cell { "library (double)", q, parallel, t });
            }
        }
    }

    return report;
}

} // namespace

TEST_CASE(
    "Barrier potential per scalar type",
    "[!benchmark][simd][simd_barrier_potential]")
{
    // The two smallest scenes have a few hundred collisions: one or two TBB
    // tasks, and timings dominated by fixed overhead rather than the kernels.
    const std::vector<std::string> skipped = { "two-cubes", "bunny" };

    // One report per available scene, in the specs' order (increasing size).
    std::vector<Report> reports;
    for (const ipc::tests::AssemblySceneSpec& spec :
         ipc::tests::assembly_scene_specs()) {
        if (std::find(skipped.begin(), skipped.end(), spec.label)
            != skipped.end()) {
            continue;
        }
        const std::optional<ipc::tests::AssemblyScene> scene =
            ipc::tests::build_assembly_scene(spec);
        if (!scene.has_value()) {
            fmt::print("Skipping {}: mesh not available\n", spec.label);
            continue;
        }
        reports.push_back(bench_scene(scene.value()));
        print_report(reports.back());
    }
    REQUIRE(!reports.empty());

    if (const char* path = std::getenv("IPC_TOOLKIT_BENCH_OUTPUT")) {
        std::ofstream f(path);
        f << "{\n  \"scenes\": [\n";
        for (size_t i = 0; i < reports.size(); i++) {
            f << reports[i].to_json()
              << (i + 1 < reports.size() ? ",\n" : "\n");
        }
        f << "  ]\n}\n";
        fmt::print("Wrote {}\n", path);
    }
}
