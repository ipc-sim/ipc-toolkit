#pragma once
#include <tests/utils.hpp>
#include <ipc/distance/distance_type.hpp>
#include <ipc/distance/point_triangle.hpp>
#include <ipc/distance/exact_vector.hpp>

using namespace ipc;
using ExReal = ExactReal;    // exact scalar type
using ExVec3 = ExactVector3; // exact vector

// Mirrors the relative sin² threshold of ipc::edge_edge_distance_type's
// parallel-edge handling (and ipc::PARALLEL_THRESHOLD).
constexpr double REFERENCE_PARALLEL_THRESHOLD = 2.5e-16;

inline void init_pck() { init_exact_arithmetic(); }

inline ExVec3 make_exact(Eigen::ConstRef<VectorMax3d> v)
{
    return ExVec3::from(v.data(), int(v.size())); // z = 0 for 2D vectors
}
inline PointEdgeDistanceType point_edge_distance_type_reference(
    Eigen::ConstRef<VectorMax3d> p_,
    Eigen::ConstRef<VectorMax3d> e0_,
    Eigen::ConstRef<VectorMax3d> e1_)
{
    init_pck();
    const ExVec3 p = make_exact(p_);
    const ExVec3 e0 = make_exact(e0_);
    const ExVec3 e1 = make_exact(e1_);
    const ExVec3 e = e1 - e0;

    if (e.dot(p - e0) <= 0) {
        return PointEdgeDistanceType::P_E0; // PP (p-e0)
    } else if ((-e).dot(p - e1) <= 0) {
        return PointEdgeDistanceType::P_E1; // PP (p-e1)
    } else {
        return PointEdgeDistanceType::P_E; // PE
    }
}

inline PointTriangleDistanceType point_triangle_distance_type_reference(
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

    if (r0.dot(e0) <= 0 && r0.dot(e2) >= 0) {
        return PointTriangleDistanceType::P_T0;
    }
    if (r1.dot(e1) <= 0 && r1.dot(e0) >= 0) {
        return PointTriangleDistanceType::P_T1;
    }
    if (r2.dot(e2) <= 0 && r2.dot(e1) >= 0) {
        return PointTriangleDistanceType::P_T2;
    }

    const ExVec3 n = e0.cross(-e1);

    if (n.dot(r0.cross(e0)) <= 0 && r0.dot(e0) > 0 && r1.dot(e0) < 0) {
        return PointTriangleDistanceType::P_E0;
    }
    if (n.dot(r1.cross(e1)) <= 0 && r1.dot(e1) > 0 && r2.dot(e1) < 0) {
        return PointTriangleDistanceType::P_E1;
    }
    if (n.dot(r2.cross(e2)) <= 0 && r2.dot(e2) > 0 && r0.dot(e2) < 0) {
        return PointTriangleDistanceType::P_E2;
    }

    return PointTriangleDistanceType::P_T;
}

inline bool is_parallel_edge_edge_reference(
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

    const ExReal cross_norm_sqr = u.cross(v).squared_norm();
    if constexpr (REFERENCE_PARALLEL_THRESHOLD == 0.0)
        return cross_norm_sqr == 0;
    const ExReal a = u.squared_norm();
    const ExReal c = v.squared_norm();
    return cross_norm_sqr < a * c * REFERENCE_PARALLEL_THRESHOLD;
}

inline EdgeEdgeDistanceType edge_edge_parallel_distance_type_reference(
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
    const ExReal ea_len_sqr = ea.squared_norm();
    const ExReal alpha_N = (eb0 - ea0).dot(ea);
    const ExReal beta_N = (eb1 - ea0).dot(ea);

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
///        edge_edge_parallel_distance_type_reference is only valid for
///        genuinely parallel edges. Defaults to the library's threshold so the
///        reference mirrors the shipped behaviour unless asked otherwise.
inline EdgeEdgeDistanceType edge_edge_distance_type_reference(
    Eigen::ConstRef<Eigen::Vector3d> ea0_,
    Eigen::ConstRef<Eigen::Vector3d> ea1_,
    Eigen::ConstRef<Eigen::Vector3d> eb0_,
    Eigen::ConstRef<Eigen::Vector3d> eb1_,
    const double parallel_threshold = REFERENCE_PARALLEL_THRESHOLD)
{
    init_pck();
    const ExVec3 ea0 = make_exact(ea0_);
    const ExVec3 ea1 = make_exact(ea1_);
    const ExVec3 eb0 = make_exact(eb0_);
    const ExVec3 eb1 = make_exact(eb1_);

    const ExVec3 u = ea1 - ea0;
    const ExVec3 v = eb1 - eb0;
    const ExVec3 w = ea0 - eb0;

    const ExReal a = u.squared_norm(); // always ≥ 0
    const ExReal b = u.dot(v);
    const ExReal c = v.squared_norm(); // always ≥ 0
    const ExReal d = u.dot(w);
    const ExReal e = v.dot(w);
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
    const ExReal cross_norm_sqr = u.cross(v).squared_norm();
    bool is_parallel;
    if (parallel_threshold == 0.0) {
        is_parallel = (cross_norm_sqr == 0);
    } else {
        is_parallel = cross_norm_sqr < a * c * parallel_threshold;
    }
    if (is_parallel) {
        return edge_edge_parallel_distance_type_reference(
            ea0_, ea1_, eb0_, eb1_);
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