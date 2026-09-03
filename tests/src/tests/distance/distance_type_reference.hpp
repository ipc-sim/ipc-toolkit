#pragma once
#include <tests/utils.hpp>
#include <ipc/distance/distance_type.hpp>
#include <ipc/distance/point_triangle.hpp>
#include <geogram/numerics/exact_geometry.h>

using namespace ipc;
using ExReal = GEO::expansion_nt; // exact scalar type
using ExVec3 = GEO::vec3E;        // exact vector

inline void init_pck()
{ // TODO init once in main
    static bool initialized = false;
    if (!initialized) {
        GEO::PCK::initialize();
        initialized = true;
    }
}
inline ExVec3 make_exact(Eigen::ConstRef<VectorMax3d> v)
{
    ExReal x { v.x() };
    ExReal y { v.y() };
    ExReal z { v.size() < 3 ? 0 : v.z() }; // compatibility with 2D vectors
    return ExVec3(std::move(x), std::move(y), std::move(z));
}
PointEdgeDistanceType point_edge_distance_type_exact(
    Eigen::ConstRef<VectorMax3d> p_,
    Eigen::ConstRef<VectorMax3d> e0_,
    Eigen::ConstRef<VectorMax3d> e1_)
{
    init_pck();
    const ExVec3 p = make_exact(p_);
    const ExVec3 e0 = make_exact(e0_);
    const ExVec3 e1 = make_exact(e1_);
    const ExVec3 e = e1 - e0;

    if (dot(e, p - e0) <= 0) {
        return PointEdgeDistanceType::P_E0; // PP (p-e0)
    } else if (dot(-e, p - e1) <= 0) {
        return PointEdgeDistanceType::P_E1; // PP (p-e1)
    } else {
        return PointEdgeDistanceType::P_E; // PE
    }
}

PointTriangleDistanceType point_triangle_distance_type_exact(
    Eigen::ConstRef<Eigen::Vector3d> p_,
    Eigen::ConstRef<Eigen::Vector3d> t0_,
    Eigen::ConstRef<Eigen::Vector3d> t1_,
    Eigen::ConstRef<Eigen::Vector3d> t2_)
{
    init_pck();
    const ExVec3 p = make_exact(p_);
    const ExVec3 t0 = make_exact(t0_);
    const ExVec3 t1 = make_exact(t1_);
    const ExVec3 t2 = make_exact(t2_);

    const ExVec3 e0 = t1 - t0;
    const ExVec3 e1 = t2 - t1;
    const ExVec3 e2 = t0 - t2;

    const ExVec3 r0 = p - t0;
    const ExVec3 r1 = p - t1;
    const ExVec3 r2 = p - t2;

    if (dot(r0, e0) <= 0 && dot(r0, e2) >= 0) {
        return PointTriangleDistanceType::P_T0;
    }
    if (dot(r1, e1) <= 0 && dot(r1, e0) >= 0) {
        return PointTriangleDistanceType::P_T1;
    }
    if (dot(r2, e2) <= 0 && dot(r2, e1) >= 0) {
        return PointTriangleDistanceType::P_T2;
    }

    const ExVec3 n = cross(e0, -e1);

    if (dot(n, cross(r0, e0)) <= 0 && dot(r0, e0) > 0 && dot(r1, e0) < 0) {
        return PointTriangleDistanceType::P_E0;
    }
    if (dot(n, cross(r1, e1)) <= 0 && dot(r1, e1) > 0 && dot(r2, e1) < 0) {
        return PointTriangleDistanceType::P_E1;
    }
    if (dot(n, cross(r2, e2)) <= 0 && dot(r2, e2) > 0 && dot(r0, e2) < 0) {
        return PointTriangleDistanceType::P_E2;
    }

    return PointTriangleDistanceType::P_T;
}

bool is_parallel_edge_edge_exact(
    Eigen::ConstRef<Eigen::Vector3d> ea0_,
    Eigen::ConstRef<Eigen::Vector3d> ea1_,
    Eigen::ConstRef<Eigen::Vector3d> eb0_,
    Eigen::ConstRef<Eigen::Vector3d> eb1_)
{
    // TODO use a zero filter?
    init_pck();
    const ExVec3 ea0 = make_exact(ea0_);
    const ExVec3 ea1 = make_exact(ea1_);
    const ExVec3 eb0 = make_exact(eb0_);
    const ExVec3 eb1 = make_exact(eb1_);

    const ExVec3 u = ea1 - ea0;
    const ExVec3 v = eb1 - eb0;

    const ExReal cross_norm_sqr = cross(u, v).length2();
    if constexpr (PARALLEL_THRESHOLD == 0.0)
        return cross_norm_sqr == 0;
    const ExReal a = u.length2();
    const ExReal c = v.length2();
    return cross_norm_sqr < a * c * PARALLEL_THRESHOLD;
}

EdgeEdgeDistanceType edge_edge_parallel_distance_type_exact(
    Eigen::ConstRef<Eigen::Vector3d> ea0_,
    Eigen::ConstRef<Eigen::Vector3d> ea1_,
    Eigen::ConstRef<Eigen::Vector3d> eb0_,
    Eigen::ConstRef<Eigen::Vector3d> eb1_)
{
    init_pck();
    const ExVec3 ea0 = make_exact(ea0_);
    const ExVec3 ea1 = make_exact(ea1_);
    const ExVec3 eb0 = make_exact(eb0_);
    const ExVec3 eb1 = make_exact(eb1_);

    const ExVec3 ea = ea1 - ea0;
    const ExReal ea_len_sqr = ea.length2();
    const ExReal alpha_N = dot(eb0 - ea0, ea);
    const ExReal beta_N = dot(eb1 - ea0, ea);

    uint8_t eac; // 0: EA0, 1: EA1, 2: EA
    uint8_t ebc; // 0: EB0, 1: EB1, 2: EB
    if (alpha_N < 0) {
        eac = (beta_N >= 0 && beta_N <= ea_len_sqr) ? 2 : 0;
        ebc = (beta_N <= alpha_N) ? 0 : (beta_N <= ea_len_sqr ? 1 : 2);
    } else if (alpha_N > ea_len_sqr) {
        eac = (beta_N >= 0 && beta_N <= ea_len_sqr) ? 2 : 1;
        ebc = (beta_N >= alpha_N) ? 0 : (beta_N >= 0 ? 1 : 2);
    } else {
        eac = 2;
        ebc = 0;
    }

    // f(0, 0) = 0000 = 0 -> EA0_EB0
    // f(0, 1) = 0001 = 1 -> EA0_EB1
    // f(1, 0) = 0010 = 2 -> EA1_EB0
    // f(1, 1) = 0011 = 3 -> EA1_EB1
    // f(2, 0) = 0100 = 4 -> EA_EB0
    // f(2, 1) = 0101 = 5 -> EA_EB1
    // f(0, 2) = 0110 = 6 -> EA0_EB
    // f(1, 2) = 0111 = 7 -> EA1_EB
    // f(2, 2) = 1000 = 8 -> EA_EB

    assert(eac != 2 || ebc != 2); // This case results in a degenerate line-line
    return EdgeEdgeDistanceType(ebc < 2 ? (eac << 1 | ebc) : (6 + eac));
}

// A more robust implementation of http://geomalgorithms.com/a07-_distance.html
/// @param parallel_threshold Relative sin² tolerance used to decide whether the
///        edges are parallel. Pass 0 for a fully exact reference: only edges
///        that are *exactly* parallel take the parallel branch. Any non-zero
///        value makes this a hybrid (exact arithmetic, approximate parallelism
///        test) which can misclassify near-parallel edges, since
///        edge_edge_parallel_distance_type_exact is only valid for genuinely
///        parallel edges. Defaults to the library's PARALLEL_THRESHOLD so the
///        reference mirrors the shipped behaviour unless asked otherwise.
EdgeEdgeDistanceType edge_edge_distance_type_exact(
    Eigen::ConstRef<Eigen::Vector3d> ea0_,
    Eigen::ConstRef<Eigen::Vector3d> ea1_,
    Eigen::ConstRef<Eigen::Vector3d> eb0_,
    Eigen::ConstRef<Eigen::Vector3d> eb1_,
    const double parallel_threshold = PARALLEL_THRESHOLD)
{
    init_pck();
    const ExVec3 ea0 = make_exact(ea0_);
    const ExVec3 ea1 = make_exact(ea1_);
    const ExVec3 eb0 = make_exact(eb0_);
    const ExVec3 eb1 = make_exact(eb1_);

    const ExVec3 u = ea1 - ea0;
    const ExVec3 v = eb1 - eb0;
    const ExVec3 w = ea0 - eb0;

    const ExReal a = u.length2(); // always ≥ 0
    const ExReal b = dot(u, v);
    const ExReal c = v.length2(); // always ≥ 0
    const ExReal d = dot(u, w);
    const ExReal e = dot(v, w);
    const ExReal D = a * c - b * b; // always ≥ 0

    // Degenerate cases should not happen in practice, but we handle them
    if (a == 0 && c == 0) {
        return EdgeEdgeDistanceType::EA0_EB0;
    } else if (a == 0) {
        return EdgeEdgeDistanceType::EA0_EB;
    } else if (c == 0) {
        return EdgeEdgeDistanceType::EA_EB0;
    }

    // Special handling for parallel edges
    const ExReal cross_norm_sqr = cross(u, v).length2();
    bool is_parallel;
    if (parallel_threshold == 0.0) {
        is_parallel = (cross_norm_sqr == 0);
    } else {
        is_parallel = cross_norm_sqr < a * c * parallel_threshold;
    }
    if (is_parallel) {
        return edge_edge_parallel_distance_type_exact(ea0_, ea1_, eb0_, eb1_);
    }

    EdgeEdgeDistanceType default_case = EdgeEdgeDistanceType::EA_EB;

    // compute the line parameters of the two closest points
    const ExReal sN = (b * e - c * d);
    ExReal tN, tD; // tc = tN / tD
    if (sN <= 0) { // sc < 0 ⟹ the s=0 edge is visible
        tN = e;
        tD = c;
        default_case = EdgeEdgeDistanceType::EA0_EB;
    } else if (sN >= D) { // sc > 1 ⟹ the s=1 edge is visible
        tN = e + b;
        tD = c;
        default_case = EdgeEdgeDistanceType::EA1_EB;
    } else {
        tN = (a * e - b * d);
        tD = D; // default tD = D ≥ 0
    }

    if (tN <= 0) { // tc < 0 ⟹ the t=0 edge is visible
        // recompute sc for this edge
        if (-d <= 0) {
            return EdgeEdgeDistanceType::EA0_EB0;
        } else if (-d >= a) {
            return EdgeEdgeDistanceType::EA1_EB0;
        } else {
            return EdgeEdgeDistanceType::EA_EB0;
        }
    } else if (tN >= tD) { // tc > 1 ⟹ the t=1 edge is visible
        // recompute sc for this edge
        if ((-d + b) <= 0) {
            return EdgeEdgeDistanceType::EA0_EB1;
        } else if ((-d + b) >= a) {
            return EdgeEdgeDistanceType::EA1_EB1;
        } else {
            return EdgeEdgeDistanceType::EA_EB1;
        }
    }

    return default_case;
}