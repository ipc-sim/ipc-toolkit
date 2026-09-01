#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>

#include <ipc/barrier/barrier.hpp>
#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/line_line.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/distance/point_line.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/distance/point_triangle.hpp>
#include <ipc/distance/edge_edge_mollifier.hpp>
#include <ipc/distance/point_plane.hpp>
#include <ipc/geometry/area.hpp>
#include <ipc/tangent/closest_point.hpp>
#include <ipc/tangent/relative_velocity.hpp>
#include <ipc/tangent/tangent_basis.hpp>

#include <Eigen/Geometry>

#include <vector>

// =============================================================================
// How the distance functions pass and bind their arguments.
//
// Every benchmark here compares the real library function ("Lib") against a
// same-TU reference implementation with the same body but a chosen parameter
// type. The reference rows double as a noise floor: they are unaffected by
// changes under src/, so if a "Lib" row moves and its control moves with it,
// the machine drifted, not the code.
//
// These are hidden ([!benchmark]); run them by name, e.g.
//     ./ipc_toolkit_tests "Parameter type (PE)"
//
// -- Method ------------------------------------------------------------------
//
// Hand-written stand-ins overstated the win three separate times during the
// v2.0.0 templating work. Only two methods proved trustworthy:
//   * stash-compare: build the same benchmark against src/ with and without
//     the change, and
//   * interleaved A/B: alternate the two binaries A/B/A/B... and take medians.
// Always report the control rows alongside the result. On a loaded machine the
// controls routinely drift 3-6%, which is larger than most effects here.
//
// Three traps this file exists to keep measured:
//
//   1. LITERAL DISTANCE TYPES. Passing EdgeEdgeDistanceType::EA_EB as a
//      literal lets the compiler fold the switch to one case after inlining,
//      which flattered the inlined version by 2-4x. The "runtime dtype" rows
//      read the type from a heap array, which is what a real caller does.
//      Measure those.
//   2. LICM. Catch2 has no DoNotOptimize, so the indices come from a
//      heap-allocated array of prime length (IDX_M): i % IDX_M has no
//      power-of-2 period, and the heap pointer stops the compiler tracing
//      provenance back to V.
//
// -- What was learned --------------------------------------------------------
//
// Final measured state of the v2.0.0 distance work, on an Apple M-series, -O3,
// call site V.row(i) of a column-major MatrixXd (the dominant in-tree
// pattern):
//
//                                         before     after
//   point_line_distance                    8.2 ns    2.3 ns   3.5x
//   point_edge_distance, AUTO             18.1 ns    5.2 ns   3.5x
//   point_triangle_distance, AUTO         93.3 ns   12.8 ns   7.3x
//   edge_edge_distance, runtime dtype      6.4 ns    6.2 ns   1.0x
//
//   1. The cost was the DYNAMIC SIZE, not the copy. Binding a Vector3d lvalue
//      by Ref or by copy is a wash (2.02 vs 2.18 ns), but a VectorMax3d
//      parameter costs 2.4x a Vector3d one. Hence the dim-templated kernels.
//   2. Never materialize into a dynamically sized intermediate. A VectorMax3d
//      local costs 5.9 ns (10.9 ns from a VectorMax3d caller). Branch on
//      size() FIRST -- which does not evaluate a lazy expression -- then copy
//      into a fixed-size local.
//   3. Moving the cold `throw` bodies out of line (detail::throw_* in
//      distance_type.hpp) was the single biggest lever, worth up to 2x by
//      itself. Constructing a std::invalid_argument inline blows the inliner's
//      budget and pushes the whole function out of line.
//   4. DO NOT INLINE A WIDE DISPATCH. Moving edge_edge_distance and
//      point_triangle_distance into their headers was tried and reverted:
//      their switches have 9 and 7 cases, and inlining them at every call site
//      loses whenever the distance type is a runtime value (6.4 -> 7.1 ns, and
//      AUTO 9.7 -> 10.5 ns). point_edge_distance has 3 cases and wins.
//   5. Eigen::ConstRef is the right default for a kernel that reads its
//      arguments more than a few times -- it guarantees a single evaluation
//      and measured 1.25-1.38x better than `const&` on point_edge's AUTO path.
//      But it is not free: binding an expression to a Ref materializes it into
//      the Ref's buffer. For line_line's gradient/Hessian, which read three
//      coefficients per argument and hand them straight to generated code,
//      that copy is pure cost -- Eigen::MatrixBase<Derived> measured 1.16x
//      faster on an expression argument and no worse anywhere else.
//   6. point_triangle_distance_type formerly cost ~87 ns, 27x the distance it
//      selects for, because it ran three 2x2 LDLT solves. It is now a closed
//      form (the Gram matrix is diagonal, since each edge is perpendicular to
//      the normal). The error study is summarized in the v2.0.0 release notes.
// =============================================================================

using namespace ipc;

namespace {

/// Random indices on the heap, to defeat LICM. Prime length, see trap 2 above.
constexpr int IDX_M = 997;
Eigen::ArrayXi random_indices(const int n)
{
    return Eigen::ArrayXi::Random(IDX_M).abs().unaryExpr(
        [n](int x) { return x % n; });
}

// --- point-point -------------------------------------------------------------

double pp_ref_3d(
    Eigen::ConstRef<Eigen::Vector3d> p0, Eigen::ConstRef<Eigen::Vector3d> p1)
{
    return (p1 - p0).squaredNorm();
}

double pp_copy_3d(const Eigen::Vector3d& p0, const Eigen::Vector3d& p1)
{
    return (p1 - p0).squaredNorm();
}

double
pp_ref_maxN(Eigen::ConstRef<VectorMax3d> p0, Eigen::ConstRef<VectorMax3d> p1)
{
    return (p1 - p0).squaredNorm();
}

// --- point-line --------------------------------------------------------------

/// @brief The shared fixed-size body. Always inlined, so any difference
/// measured below comes from the parameter type, not from this.
template <int N, typename DerivedP, typename DerivedE0, typename DerivedE1>
inline double
pl_body(const DerivedP& p, const DerivedE0& e0, const DerivedE1& e1)
{
    if constexpr (N == 2) {
        const Eigen::Vector2d e = e1 - e0;
        const double numerator =
            e[1] * p[0] - e[0] * p[1] + e1[0] * e0[1] - e1[1] * e0[0];
        return numerator * numerator / e.squaredNorm();
    } else {
        const Eigen::Vector3d p_to_e0 = e0 - p, p_to_e1 = e1 - p;
        return p_to_e0.cross(p_to_e1).squaredNorm() / (e1 - e0).squaredNorm();
    }
}

/// The pre-v2.0.0 signature for the dimension-agnostic functions: dynamic Ref.
double pl_ref_maxN(
    Eigen::ConstRef<VectorMax3d> p,
    Eigen::ConstRef<VectorMax3d> e0,
    Eigen::ConstRef<VectorMax3d> e1)
{
    if (p.size() == 2) {
        return pl_body<2>(p.head<2>(), e0.head<2>(), e1.head<2>());
    }
    return pl_body<3>(p.head<3>(), e0.head<3>(), e1.head<3>());
}

/// The best this body can do: a statically sized Ref. The bar for "Lib".
double pl_ref_3d(
    Eigen::ConstRef<Eigen::Vector3d> p,
    Eigen::ConstRef<Eigen::Vector3d> e0,
    Eigen::ConstRef<Eigen::Vector3d> e1)
{
    return pl_body<3>(p, e0, e1);
}

// --- point-edge: the same question, with the distance-type dispatch on top ---

/// @brief Fixed-size stand-in for point_edge_distance_type. Matches the
/// library apart from the degenerate-edge logger call, which never fires here.
template <int N, typename DerivedP, typename DerivedE0, typename DerivedE1>
inline PointEdgeDistanceType
pe_type_body(const DerivedP& p, const DerivedE0& e0, const DerivedE1& e1)
{
    const Eigen::Vector<double, N> e = e1 - e0;
    const double e_length_sqr = e.squaredNorm();
    if (e_length_sqr == 0) {
        return PointEdgeDistanceType::P_E0;
    }
    const double ratio = e.dot(p - e0) / e_length_sqr;
    if (ratio < 0) {
        return PointEdgeDistanceType::P_E0;
    } else if (ratio > 1) {
        return PointEdgeDistanceType::P_E1;
    }
    return PointEdgeDistanceType::P_E;
}

/// @brief The shared fixed-size body, as pl_body but with the dtype dispatch.
template <int N, typename DerivedP, typename DerivedE0, typename DerivedE1>
inline double pe_body(
    const DerivedP& p,
    const DerivedE0& e0,
    const DerivedE1& e1,
    PointEdgeDistanceType dtype)
{
    if (dtype == PointEdgeDistanceType::AUTO) {
        dtype = pe_type_body<N>(p, e0, e1);
    }
    switch (dtype) {
    case PointEdgeDistanceType::P_E0:
        return (p - e0).squaredNorm();
    case PointEdgeDistanceType::P_E1:
        return (p - e1).squaredNorm();
    default:
        return pl_body<N>(p, e0, e1);
    }
}

/// The pre-v2.0.0 signature: a dynamic-size Ref.
double pe_ref_maxN(
    Eigen::ConstRef<VectorMax3d> p,
    Eigen::ConstRef<VectorMax3d> e0,
    Eigen::ConstRef<VectorMax3d> e1,
    PointEdgeDistanceType dtype)
{
    if (p.size() == 2) {
        return pe_body<2>(p.head<2>(), e0.head<2>(), e1.head<2>(), dtype);
    }
    return pe_body<3>(p.head<3>(), e0.head<3>(), e1.head<3>(), dtype);
}

/// The bar for "Lib": a statically sized Ref.
double pe_ref_3d(
    Eigen::ConstRef<Eigen::Vector3d> p,
    Eigen::ConstRef<Eigen::Vector3d> e0,
    Eigen::ConstRef<Eigen::Vector3d> e1,
    PointEdgeDistanceType dtype)
{
    return pe_body<3>(p, e0, e1, dtype);
}

// --- line-line ---------------------------------------------------------------

double ll_ref_3d(
    Eigen::ConstRef<Eigen::Vector3d> ea0,
    Eigen::ConstRef<Eigen::Vector3d> ea1,
    Eigen::ConstRef<Eigen::Vector3d> eb0,
    Eigen::ConstRef<Eigen::Vector3d> eb1)
{
    const Eigen::Vector3d normal = (ea1 - ea0).cross(eb1 - eb0);
    const double line_to_line = (eb0 - ea0).dot(normal);
    return line_to_line * line_to_line / normal.squaredNorm();
}

/// A Ref that can bind a row of a column-major matrix WITHOUT copying, by
/// admitting a runtime inner stride. The plain Ref<const Vector3d> above
/// cannot, so it materializes the row into its own buffer first.
double ll_row_ref_3d(
    const Eigen::Ref<const Eigen::Vector3d, 0, Eigen::InnerStride<>>& ea0,
    const Eigen::Ref<const Eigen::Vector3d, 0, Eigen::InnerStride<>>& ea1,
    const Eigen::Ref<const Eigen::Vector3d, 0, Eigen::InnerStride<>>& eb0,
    const Eigen::Ref<const Eigen::Vector3d, 0, Eigen::InnerStride<>>& eb1)
{
    const Eigen::Vector3d normal = (ea1 - ea0).cross(eb1 - eb0);
    const double line_to_line = (eb0 - ea0).dot(normal);
    return line_to_line * line_to_line / normal.squaredNorm();
}

double ll_copy_3d(
    const Eigen::Vector3d& ea0,
    const Eigen::Vector3d& ea1,
    const Eigen::Vector3d& eb0,
    const Eigen::Vector3d& eb1)
{
    const Eigen::Vector3d normal = (ea1 - ea0).cross(eb1 - eb0);
    const double line_to_line = (eb0 - ea0).dot(normal);
    return line_to_line * line_to_line / normal.squaredNorm();
}

Eigen::Vector<double, 12> ll_grad_ref_3d(
    Eigen::ConstRef<Eigen::Vector3d> ea0,
    Eigen::ConstRef<Eigen::Vector3d> ea1,
    Eigen::ConstRef<Eigen::Vector3d> eb0,
    Eigen::ConstRef<Eigen::Vector3d> eb1)
{
    Eigen::Vector<double, 12> grad;
    autogen::line_line_distance_gradient(
        ea0[0], ea0[1], ea0[2], ea1[0], ea1[1], ea1[2], eb0[0], eb0[1], eb0[2],
        eb1[0], eb1[1], eb1[2], grad.data());
    return grad;
}

Eigen::Matrix<double, 12, 12> ll_hess_ref_3d(
    Eigen::ConstRef<Eigen::Vector3d> ea0,
    Eigen::ConstRef<Eigen::Vector3d> ea1,
    Eigen::ConstRef<Eigen::Vector3d> eb0,
    Eigen::ConstRef<Eigen::Vector3d> eb1)
{
    Eigen::Matrix<double, 12, 12> hess;
    autogen::line_line_distance_hessian(
        ea0[0], ea0[1], ea0[2], ea1[0], ea1[1], ea1[2], eb0[0], eb0[1], eb0[2],
        eb1[0], eb1[1], eb1[2], hess.data());
    return hess;
}

/// Writes through an out-parameter instead of returning. Paired with the
/// returning form above to check that NRVO makes the two equivalent.
template <typename DerivedHess>
void ll_hess_out_param(
    Eigen::ConstRef<Eigen::Vector3d> ea0,
    Eigen::ConstRef<Eigen::Vector3d> ea1,
    Eigen::ConstRef<Eigen::Vector3d> eb0,
    Eigen::ConstRef<Eigen::Vector3d> eb1,
    Eigen::PlainObjectBase<DerivedHess>& hess)
{
    hess.resize(12, 12);
    autogen::line_line_distance_hessian(
        ea0[0], ea0[1], ea0[2], ea1[0], ea1[1], ea1[2], eb0[0], eb0[1], eb0[2],
        eb1[0], eb1[1], eb1[2], hess.data());
}

// --- edge-edge ---------------------------------------------------------------

/// A same-TU, statically sized reference implementation: the best this body
/// can do, against which the library function is judged.
double ee_ref_3d(
    Eigen::ConstRef<Eigen::Vector3d> ea0,
    Eigen::ConstRef<Eigen::Vector3d> ea1,
    Eigen::ConstRef<Eigen::Vector3d> eb0,
    Eigen::ConstRef<Eigen::Vector3d> eb1,
    EdgeEdgeDistanceType dtype)
{
    if (dtype == EdgeEdgeDistanceType::AUTO) {
        dtype = edge_edge_distance_type(ea0, ea1, eb0, eb1);
    }
    switch (dtype) {
    case EdgeEdgeDistanceType::EA0_EB0:
        return (ea0 - eb0).squaredNorm();
    case EdgeEdgeDistanceType::EA0_EB1:
        return (ea0 - eb1).squaredNorm();
    case EdgeEdgeDistanceType::EA1_EB0:
        return (ea1 - eb0).squaredNorm();
    case EdgeEdgeDistanceType::EA1_EB1:
        return (ea1 - eb1).squaredNorm();
    case EdgeEdgeDistanceType::EA_EB0:
        return pl_body<3>(eb0, ea0, ea1);
    case EdgeEdgeDistanceType::EA_EB1:
        return pl_body<3>(eb1, ea0, ea1);
    case EdgeEdgeDistanceType::EA0_EB:
        return pl_body<3>(ea0, eb0, eb1);
    case EdgeEdgeDistanceType::EA1_EB:
        return pl_body<3>(ea1, eb0, eb1);
    default: {
        const Eigen::Vector3d n = (ea1 - ea0).cross(eb1 - eb0);
        const double ltl = (eb0 - ea0).dot(n);
        return ltl * ltl / n.squaredNorm();
    }
    }
}

} // namespace

// Rows of a column-major matrix: the dominant in-tree call site. The _D
// variants take a trailing distance type. (No __VA_OPT__ -- this is C++17.)
#define ROWS3(f)                                                               \
    f(V.row(idx[i % IDX_M]), V.row(idx[(i + 1) % IDX_M]),                      \
      V.row(idx[(i + 2) % IDX_M]))
#define ROWS3_D(f, dt)                                                         \
    f(V.row(idx[i % IDX_M]), V.row(idx[(i + 1) % IDX_M]),                      \
      V.row(idx[(i + 2) % IDX_M]), dt)
#define ROWS4(f)                                                               \
    f(V.row(idx[i % IDX_M]), V.row(idx[(i + 1) % IDX_M]),                      \
      V.row(idx[(i + 2) % IDX_M]), V.row(idx[(i + 3) % IDX_M]))
#define ROWS4_D(f, dt)                                                         \
    f(V.row(idx[i % IDX_M]), V.row(idx[(i + 1) % IDX_M]),                      \
      V.row(idx[(i + 2) % IDX_M]), V.row(idx[(i + 3) % IDX_M]), dt)

// The same call site, but from storage whose type already knows the dimension.
#define AT3(f, src)                                                            \
    f(src[idx[i % IDX_M]], src[idx[(i + 1) % IDX_M]], src[idx[(i + 2) % IDX_M]])
#define AT3_D(f, src, dt)                                                      \
    f(src[idx[i % IDX_M]], src[idx[(i + 1) % IDX_M]],                          \
      src[idx[(i + 2) % IDX_M]], dt)
#define AT4(f, src)                                                            \
    f(src[idx[i % IDX_M]], src[idx[(i + 1) % IDX_M]],                          \
      src[idx[(i + 2) % IDX_M]], src[idx[(i + 3) % IDX_M]])
#define AT4_D(f, src, dt)                                                      \
    f(src[idx[i % IDX_M]], src[idx[(i + 1) % IDX_M]],                          \
      src[idx[(i + 2) % IDX_M]], src[idx[(i + 3) % IDX_M]], dt)

TEST_CASE("Parameter type (PP)", "[!benchmark][eigen][pp][param]")
{
    constexpr int N = 100;
    const Eigen::MatrixXd V = Eigen::MatrixXd::Random(N, 3);
    const Eigen::ArrayXi idx = random_indices(N);

#define PP_ROWS(f) f(V.row(idx[i % IDX_M]), V.row(idx[(i + 1) % IDX_M]))

    BENCHMARK("rows: Lib", i) { return PP_ROWS(point_point_distance); };
    BENCHMARK("rows: Ref<Vector3d>", i) { return PP_ROWS(pp_ref_3d); };
    BENCHMARK("rows: Ref<VectorMax3d>", i) { return PP_ROWS(pp_ref_maxN); };
    BENCHMARK("rows: by-value Vector3d", i) { return PP_ROWS(pp_copy_3d); };

#undef PP_ROWS
}

TEST_CASE("Parameter type (PL)", "[!benchmark][eigen][pl][param]")
{
    constexpr int N = 100;
    const Eigen::MatrixXd V = Eigen::MatrixXd::Random(N, 3);
    const Eigen::ArrayXi idx = random_indices(N);

    // Call sites that already know the dimension statically, and ones that
    // reach the front end's runtime branch.
    std::vector<Eigen::Vector3d> P3(N);
    std::vector<VectorMax3d> PN(N);
    for (int k = 0; k < N; ++k) {
        P3[k] = V.row(k);
        PN[k] = V.row(k);
    }

    BENCHMARK("rows: Lib", i) { return ROWS3(point_line_distance); };
    BENCHMARK("rows: Ref<Vector3d>", i) { return ROWS3(pl_ref_3d); };
    BENCHMARK("rows: Ref<VectorMax3d>", i) { return ROWS3(pl_ref_maxN); };

    BENCHMARK("Vector3d: Lib", i) { return AT3(point_line_distance, P3); };
    BENCHMARK("Vector3d: Ref<Vector3d>", i) { return AT3(pl_ref_3d, P3); };
    BENCHMARK("VectorMax3d: Lib", i) { return AT3(point_line_distance, PN); };
    BENCHMARK("VectorMax3d: Ref<VectorMax3d>", i)
    {
        return AT3(pl_ref_maxN, PN);
    };
}

TEST_CASE("Parameter type (PE)", "[!benchmark][eigen][pe][param]")
{
    constexpr int N = 100;
    const Eigen::MatrixXd V = Eigen::MatrixXd::Random(N, 3);
    const Eigen::ArrayXi idx = random_indices(N);

    std::vector<Eigen::Vector3d> P3(N);
    std::vector<VectorMax3d> PN(N);
    for (int k = 0; k < N; ++k) {
        P3[k] = V.row(k);
        PN[k] = V.row(k);
    }

    // A runtime-varying dtype: a real caller gets this from a stencil, so the
    // switch cannot fold to one case the way a literal lets it. See trap 1.
    std::vector<PointEdgeDistanceType> dts(IDX_M);
    for (int k = 0; k < IDX_M; ++k) {
        dts[k] = static_cast<PointEdgeDistanceType>(k % 3);
    }

    // AUTO: the distance type is resolved inside, as in normal use. This pays
    // the parameter-type penalty twice -- once in point_edge_distance_type and
    // once in the distance itself -- so it is the strongest signal here.
    BENCHMARK("AUTO rows: Lib", i)
    {
        return ROWS3_D(point_edge_distance, PointEdgeDistanceType::AUTO);
    };
    BENCHMARK("AUTO rows: Ref<Vector3d>", i)
    {
        return ROWS3_D(pe_ref_3d, PointEdgeDistanceType::AUTO);
    };
    BENCHMARK("AUTO rows: Ref<VectorMax3d>", i)
    {
        return ROWS3_D(pe_ref_maxN, PointEdgeDistanceType::AUTO);
    };

    BENCHMARK("runtime dtype rows: Lib", i)
    {
        return ROWS3_D(point_edge_distance, dts[i % IDX_M]);
    };
    BENCHMARK("runtime dtype rows: Ref<Vector3d>", i)
    {
        return ROWS3_D(pe_ref_3d, dts[i % IDX_M]);
    };

    BENCHMARK("AUTO Vector3d: Lib", i)
    {
        return AT3_D(point_edge_distance, P3, PointEdgeDistanceType::AUTO);
    };
    BENCHMARK("AUTO VectorMax3d: Lib", i)
    {
        return AT3_D(point_edge_distance, PN, PointEdgeDistanceType::AUTO);
    };
    BENCHMARK("AUTO VectorMax3d: Ref<VectorMax3d>", i)
    {
        return AT3_D(pe_ref_maxN, PN, PointEdgeDistanceType::AUTO);
    };

    // An argument the front end must not evaluate more than once.
    BENCHMARK("AUTO lazy expr: Lib", i)
    {
        return point_edge_distance(
            V.row(idx[i % IDX_M]) - V.row(idx[(i + 3) % IDX_M]),
            V.row(idx[(i + 1) % IDX_M]), V.row(idx[(i + 2) % IDX_M]),
            PointEdgeDistanceType::AUTO);
    };
}

TEST_CASE("Parameter type (LL)", "[!benchmark][eigen][ll][param]")
{
    constexpr int N = 100;
    const Eigen::MatrixXd V = Eigen::MatrixXd::Random(N, 3);
    const Eigen::Matrix<double, -1, -1, Eigen::RowMajor> VR = V;
    const Eigen::ArrayXi idx = random_indices(N);

    std::vector<Eigen::Vector3d> P3(N);
    for (int k = 0; k < N; ++k) {
        P3[k] = V.row(k);
    }

    BENCHMARK("rows: Lib", i) { return ROWS4(line_line_distance); };
    BENCHMARK("rows: Ref<Vector3d>", i) { return ROWS4(ll_ref_3d); };
    BENCHMARK("rows: Ref<Vector3d, InnerStride>", i)
    {
        return ROWS4(ll_row_ref_3d);
    };
    BENCHMARK("rows: by-value Vector3d", i) { return ROWS4(ll_copy_3d); };

    // A row-major matrix makes each row contiguous, so Ref<const Vector3d> can
    // bind it without materializing.
    BENCHMARK("row-major rows: Lib", i)
    {
        return line_line_distance(
            VR.row(idx[i % IDX_M]), VR.row(idx[(i + 1) % IDX_M]),
            VR.row(idx[(i + 2) % IDX_M]), VR.row(idx[(i + 3) % IDX_M]));
    };

    BENCHMARK("Vector3d: Lib", i) { return AT4(line_line_distance, P3); };
    BENCHMARK("Vector3d: Ref<Vector3d>", i) { return AT4(ll_ref_3d, P3); };
}

TEST_CASE("Parameter type (EE)", "[!benchmark][eigen][ee][param]")
{
    constexpr int N = 100;
    const Eigen::MatrixXd V = Eigen::MatrixXd::Random(N, 3);
    const Eigen::ArrayXi idx = random_indices(N);

    std::vector<Eigen::Vector3d> P3(N);
    for (int k = 0; k < N; ++k) {
        P3[k] = V.row(k);
    }

    // See trap 1: a literal dtype folds the switch, a runtime one does not.
    std::vector<EdgeEdgeDistanceType> dts(IDX_M);
    for (int k = 0; k < IDX_M; ++k) {
        dts[k] = static_cast<EdgeEdgeDistanceType>(k % 9);
    }

    BENCHMARK("AUTO rows: Lib", i)
    {
        return ROWS4_D(edge_edge_distance, EdgeEdgeDistanceType::AUTO);
    };
    BENCHMARK("AUTO rows: Ref<Vector3d>", i)
    {
        return ROWS4_D(ee_ref_3d, EdgeEdgeDistanceType::AUTO);
    };
    BENCHMARK("AUTO Vector3d: Lib", i)
    {
        return AT4_D(edge_edge_distance, P3, EdgeEdgeDistanceType::AUTO);
    };
    BENCHMARK("AUTO Vector3d: Ref<Vector3d>", i)
    {
        return AT4_D(ee_ref_3d, P3, EdgeEdgeDistanceType::AUTO);
    };

    BENCHMARK("runtime dtype rows: Lib", i)
    {
        return ROWS4_D(edge_edge_distance, dts[i % IDX_M]);
    };
    BENCHMARK("runtime dtype rows: Ref<Vector3d>", i)
    {
        return ROWS4_D(ee_ref_3d, dts[i % IDX_M]);
    };

    // Explicit dtypes isolate the parameter/inlining question from the
    // distance-type resolution. EA_EB is the interior-interior case.
    BENCHMARK("EA_EB rows: Lib", i)
    {
        return ROWS4_D(edge_edge_distance, EdgeEdgeDistanceType::EA_EB);
    };
    BENCHMARK("EA_EB rows: Ref<Vector3d>", i)
    {
        return ROWS4_D(ee_ref_3d, EdgeEdgeDistanceType::EA_EB);
    };
    BENCHMARK("EA0_EB0 rows: Lib", i)
    {
        return ROWS4_D(edge_edge_distance, EdgeEdgeDistanceType::EA0_EB0);
    };
    BENCHMARK("EA_EB0 rows: Lib", i)
    {
        return ROWS4_D(edge_edge_distance, EdgeEdgeDistanceType::EA_EB0);
    };
}

TEST_CASE("Parameter type (PT)", "[!benchmark][eigen][pt][param]")
{
    constexpr int N = 100;
    const Eigen::MatrixXd V = Eigen::MatrixXd::Random(N, 3);
    const Eigen::ArrayXi idx = random_indices(N);

    std::vector<Eigen::Vector3d> P3(N);
    for (int k = 0; k < N; ++k) {
        P3[k] = V.row(k);
    }

    BENCHMARK("AUTO rows: Lib", i)
    {
        return ROWS4_D(
            point_triangle_distance, PointTriangleDistanceType::AUTO);
    };
    BENCHMARK("P_T rows: Lib", i)
    {
        return ROWS4_D(point_triangle_distance, PointTriangleDistanceType::P_T);
    };
    BENCHMARK("P_E0 rows: Lib", i)
    {
        return ROWS4_D(
            point_triangle_distance, PointTriangleDistanceType::P_E0);
    };
    BENCHMARK("AUTO Vector3d: Lib", i)
    {
        return AT4_D(
            point_triangle_distance, P3, PointTriangleDistanceType::AUTO);
    };

    // Static vs dynamic argument type on the Hessian, which is the heaviest
    // consumer of the fixed-size return types.
    BENCHMARK("hessian rows: Lib", i)
    {
        return ROWS4(point_triangle_distance_hessian);
    };
    BENCHMARK("hessian Vector3d: Lib", i)
    {
        return AT4(point_triangle_distance_hessian, P3);
    };
}

TEST_CASE("Parameter type (LL derivatives)", "[!benchmark][eigen][ll][deriv]")
{
    constexpr int N = 100;
    const Eigen::MatrixXd V = Eigen::MatrixXd::Random(N, 3);
    const Eigen::ArrayXi idx = random_indices(N);

    std::vector<Eigen::Vector3d> P3(N);
    for (int k = 0; k < N; ++k) {
        P3[k] = V.row(k);
    }

    BENCHMARK("grad rows: Lib", i)
    {
        return ROWS4(line_line_distance_gradient);
    };
    BENCHMARK("grad rows: Ref<Vector3d>", i) { return ROWS4(ll_grad_ref_3d); };
    BENCHMARK("grad Vector3d: Lib", i)
    {
        return AT4(line_line_distance_gradient, P3);
    };
    BENCHMARK("grad Vector3d: Ref<Vector3d>", i)
    {
        return AT4(ll_grad_ref_3d, P3);
    };

    BENCHMARK("hess rows: Lib", i)
    {
        return ROWS4(line_line_distance_hessian);
    };
    BENCHMARK("hess rows: Ref<Vector3d>", i) { return ROWS4(ll_hess_ref_3d); };
    BENCHMARK("hess Vector3d: Lib", i)
    {
        return AT4(line_line_distance_hessian, P3);
    };
    BENCHMARK("hess Vector3d: Ref<Vector3d>", i)
    {
        return AT4(ll_hess_ref_3d, P3);
    };
    BENCHMARK("hess Vector3d: out-parameter", i)
    {
        Matrix12d hess;
        ll_hess_out_param(
            P3[idx[i % IDX_M]], P3[idx[(i + 1) % IDX_M]],
            P3[idx[(i + 2) % IDX_M]], P3[idx[(i + 3) % IDX_M]], hess);
        return hess;
    };

    // The real in-library caller: edge_edge_distance_{gradient,hessian} routes
    // the EA_EB (interior-interior) case straight to line_line.
    BENCHMARK("EE grad EA_EB rows: Lib", i)
    {
        return ROWS4_D(
            edge_edge_distance_gradient, EdgeEdgeDistanceType::EA_EB);
    };
    BENCHMARK("EE hess EA_EB rows: Lib", i)
    {
        return ROWS4_D(edge_edge_distance_hessian, EdgeEdgeDistanceType::EA_EB);
    };

    // These read three coefficients per argument, so a non-trivial expression
    // is re-evaluated three times rather than materialized once. Measured: a
    // cheap elementwise expression costs nothing, and even a 3x3 product is
    // within noise of materializing first -- which is why these two take
    // Eigen::MatrixBase<Derived> rather than Eigen::ConstRef. See finding 5.
    const Eigen::Matrix3d M = Eigen::Matrix3d::Random();

    BENCHMARK("grad expr diff: Lib", i)
    {
        return line_line_distance_gradient(
            V.row(idx[i % IDX_M]) - V.row(idx[(i + 4) % IDX_M]),
            V.row(idx[(i + 1) % IDX_M]), V.row(idx[(i + 2) % IDX_M]),
            V.row(idx[(i + 3) % IDX_M]));
    };
    BENCHMARK("grad expr diff: materialized", i)
    {
        const Eigen::Vector3d d =
            V.row(idx[i % IDX_M]) - V.row(idx[(i + 4) % IDX_M]);
        return line_line_distance_gradient(
            d, V.row(idx[(i + 1) % IDX_M]), V.row(idx[(i + 2) % IDX_M]),
            V.row(idx[(i + 3) % IDX_M]));
    };
    BENCHMARK("grad expr product: Lib", i)
    {
        return line_line_distance_gradient(
            M * P3[idx[i % IDX_M]], P3[idx[(i + 1) % IDX_M]],
            P3[idx[(i + 2) % IDX_M]], P3[idx[(i + 3) % IDX_M]]);
    };
    BENCHMARK("grad expr product: materialized", i)
    {
        const Eigen::Vector3d d = M * P3[idx[i % IDX_M]];
        return line_line_distance_gradient(
            d, P3[idx[(i + 1) % IDX_M]], P3[idx[(i + 2) % IDX_M]],
            P3[idx[(i + 3) % IDX_M]]);
    };
}

TEST_CASE("Float vs double", "[!benchmark][eigen][float]")
{
    constexpr int N = 100;
    const Eigen::MatrixXf Vf = Eigen::MatrixXf::Random(N, 3);
    const Eigen::MatrixXd Vd = Vf.cast<double>();
    const Eigen::ArrayXi idx = random_indices(N);

#define EE_ROWS_OF(Mat)                                                        \
    Mat.row(idx[i % IDX_M]), Mat.row(idx[(i + 1) % IDX_M]),                    \
        Mat.row(idx[(i + 2) % IDX_M]), Mat.row(idx[(i + 3) % IDX_M])

    BENCHMARK("distance: float", i)
    {
        return barrier(edge_edge_distance(EE_ROWS_OF(Vf)), 1.0f);
    };
    BENCHMARK("distance: double", i)
    {
        return barrier(edge_edge_distance(EE_ROWS_OF(Vd)), 1.0);
    };

    // The full barrier Hessian assembly, which is what a simulator actually
    // spends its time on.
    BENCHMARK("barrier hessian: float", i)
    {
        const auto d = edge_edge_distance(EE_ROWS_OF(Vf));
        const auto grad_d = edge_edge_distance_gradient(EE_ROWS_OF(Vf));
        const auto hess_d = edge_edge_distance_hessian(EE_ROWS_OF(Vf));
        const float db = barrier_first_derivative<float>(d, 1.0f);
        const float d2b = barrier_second_derivative<float>(d, 1.0f);
        return (db * hess_d + (d2b * grad_d) * grad_d.transpose()).eval();
    };
    BENCHMARK("barrier hessian: double", i)
    {
        const auto d = edge_edge_distance(EE_ROWS_OF(Vd));
        const auto grad_d = edge_edge_distance_gradient(EE_ROWS_OF(Vd));
        const auto hess_d = edge_edge_distance_hessian(EE_ROWS_OF(Vd));
        const double db = barrier_first_derivative<double>(d, 1.0);
        const double d2b = barrier_second_derivative<double>(d, 1.0);
        return (db * hess_d + (d2b * grad_d) * grad_d.transpose()).eval();
    };

#undef EE_ROWS_OF
}

#undef ROWS3
#undef ROWS3_D
#undef ROWS4
#undef ROWS4_D
#undef AT3
#undef AT3_D
#undef AT4
#undef AT4_D

// =============================================================================
// Converted families: tangent basis, closest point, relative velocity, area,
// edge-edge mollifier, point-plane. Measured old (dim-erased, double-only)
// vs. new (two-layer, dim-templated) via interleaved A/B of two binaries.
// The same-TU rows are identical in both binaries: the noise-floor control.

namespace {

double pe_cp_ref_3d(
    Eigen::ConstRef<Eigen::Vector3d> p,
    Eigen::ConstRef<Eigen::Vector3d> e0,
    Eigen::ConstRef<Eigen::Vector3d> e1)
{
    const Eigen::Vector3d e = e1 - e0;
    return e.dot(p - e0) / e.squaredNorm();
}

double el_ref_3d(
    Eigen::ConstRef<Eigen::Vector3d> e0, Eigen::ConstRef<Eigen::Vector3d> e1)
{
    return (e1 - e0).norm();
}

} // namespace

TEST_CASE("Converted families", "[!benchmark][eigen][converted]")
{
    constexpr int N = 100;
    const Eigen::MatrixXd V = Eigen::MatrixXd::Random(N, 3);
    const Eigen::ArrayXi idx = random_indices(N);

    std::vector<Eigen::Vector3d> P3(N);
    std::vector<VectorMax3d> PN(N);
    for (int k = 0; k < N; ++k) {
        P3[k] = V.row(k);
        PN[k] = V.row(k);
    }

#define R2(f) f(V.row(idx[i % IDX_M]), V.row(idx[(i + 1) % IDX_M]))
#define R3(f)                                                                  \
    f(V.row(idx[i % IDX_M]), V.row(idx[(i + 1) % IDX_M]),                      \
      V.row(idx[(i + 2) % IDX_M]))
#define R4(f)                                                                  \
    f(V.row(idx[i % IDX_M]), V.row(idx[(i + 1) % IDX_M]),                      \
      V.row(idx[(i + 2) % IDX_M]), V.row(idx[(i + 3) % IDX_M]))
#define M3(f, src)                                                             \
    f(src[idx[i % IDX_M]], src[idx[(i + 1) % IDX_M]], src[idx[(i + 2) % IDX_M]])

    // --- controls (same code in both binaries) ---
    BENCHMARK("ctl pe_closest_point Ref<Vector3d>", i)
    {
        return R3(pe_cp_ref_3d);
    };
    BENCHMARK("ctl edge_length Ref<Vector3d>", i) { return R2(el_ref_3d); };

    // --- closest point ---
    BENCHMARK("pe_closest_point rows", i)
    {
        return R3(point_edge_closest_point);
    };
    BENCHMARK("pe_closest_point Vector3d", i)
    {
        return M3(point_edge_closest_point, P3);
    };
    BENCHMARK("pe_closest_point VectorMax3d", i)
    {
        return M3(point_edge_closest_point, PN);
    };
    BENCHMARK("pe_closest_point_jacobian rows", i)
    {
        return R3(point_edge_closest_point_jacobian);
    };
    BENCHMARK("ee_closest_point rows", i)
    {
        return R4(edge_edge_closest_point);
    };
    BENCHMARK("pt_closest_point rows", i)
    {
        return R4(point_triangle_closest_point);
    };

    // --- tangent basis ---
    BENCHMARK("pp_tangent_basis rows", i)
    {
        return R2(point_point_tangent_basis);
    };
    BENCHMARK("pp_tangent_basis Vector3d", i)
    {
        return point_point_tangent_basis(
            P3[idx[i % IDX_M]], P3[idx[(i + 1) % IDX_M]]);
    };
    BENCHMARK("pe_tangent_basis rows", i)
    {
        return R3(point_edge_tangent_basis);
    };
    BENCHMARK("ee_tangent_basis rows", i)
    {
        return R4(edge_edge_tangent_basis);
    };
    BENCHMARK("pt_tangent_basis rows", i)
    {
        return R4(point_triangle_tangent_basis);
    };

    // --- relative velocity ---
    BENCHMARK("pp_relative_velocity rows", i)
    {
        return R2(point_point_relative_velocity);
    };
    BENCHMARK("pp_relative_velocity VectorMax3d", i)
    {
        return point_point_relative_velocity(
            PN[idx[i % IDX_M]], PN[idx[(i + 1) % IDX_M]]);
    };
    BENCHMARK("pp_rel_vel_jacobian(dim=3)", i)
    {
        return point_point_relative_velocity_jacobian(2 + (idx[i % IDX_M] & 1));
    };

    // --- area ---
    BENCHMARK("edge_length rows", i) { return R2(edge_length); };
    BENCHMARK("edge_length VectorMax3d", i)
    {
        return edge_length(PN[idx[i % IDX_M]], PN[idx[(i + 1) % IDX_M]]);
    };
    BENCHMARK("triangle_area rows", i) { return R3(triangle_area); };

    // --- mollifier / point-plane ---
    BENCHMARK("ee_cross_squarednorm rows", i)
    {
        return R4(edge_edge_cross_squarednorm);
    };
    BENCHMARK("point_plane_distance rows", i)
    {
        return R3(point_plane_distance);
    };

#undef R2
#undef R3
#undef R4
#undef M3
}
