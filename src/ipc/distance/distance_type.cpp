#include "distance_type.hpp"

#include <ipc/utils/eigen_ext.hpp>
#include <ipc/utils/logger.hpp>

#include <Eigen/Geometry>

#include <limits>

namespace ipc {

template <typename T>
PointEdgeDistanceType point_edge_distance_type(
    Eigen::ConstRef<VectorMax3<T>> p,
    Eigen::ConstRef<VectorMax3<T>> e0,
    Eigen::ConstRef<VectorMax3<T>> e1)
{
    assert(p.size() == 2 || p.size() == 3);
    assert(e0.size() == 2 || e0.size() == 3);
    assert(e1.size() == 2 || e1.size() == 3);

    const VectorMax3<T> e = e1 - e0;
    const T e_length_sqr = e.squaredNorm();
    if (e_length_sqr == 0) {
        logger().warn("Degenerate edge in point_edge_distance_type!");
        return PointEdgeDistanceType::P_E0; // WARNING: use arbitrary end-point
    }

    const T ratio = e.dot(p - e0) / e_length_sqr;
    if (ratio < 0) {
        return PointEdgeDistanceType::P_E0; // PP (p-e0)
    } else if (ratio > 1) {
        return PointEdgeDistanceType::P_E1; // PP (p-e1)
    } else {
        return PointEdgeDistanceType::P_E; // PE
    }
}

template <typename T>
PointTriangleDistanceType point_triangle_distance_type(
    Eigen::ConstRef<Eigen::Vector3<T>> p,
    Eigen::ConstRef<Eigen::Vector3<T>> t0,
    Eigen::ConstRef<Eigen::Vector3<T>> t1,
    Eigen::ConstRef<Eigen::Vector3<T>> t2)
{
    const Eigen::Vector3<T> normal = (t1 - t0).cross(t2 - t0);

    Eigen::Matrix<T, 2, 3> basis, param;

    basis.row(0) = t1 - t0;
    basis.row(1) = basis.row(0).cross(normal);
    param.col(0) = (basis * basis.transpose()).ldlt().solve(basis * (p - t0));
    if (param(0, 0) > 0 && param(0, 0) < 1 && param(1, 0) >= 0) {
        return PointTriangleDistanceType::P_E0; // edge 0 is the closest
    }

    basis.row(0) = t2 - t1;
    basis.row(1) = basis.row(0).cross(normal);
    param.col(1) = (basis * basis.transpose()).ldlt().solve(basis * (p - t1));
    if (param(0, 1) > 0 && param(0, 1) < 1 && param(1, 1) >= 0) {
        return PointTriangleDistanceType::P_E1; // edge 1 is the closest
    }

    basis.row(0) = t0 - t2;
    basis.row(1) = basis.row(0).cross(normal);
    param.col(2) = (basis * basis.transpose()).ldlt().solve(basis * (p - t2));
    if (param(0, 2) > 0 && param(0, 2) < 1 && param(1, 2) >= 0) {
        return PointTriangleDistanceType::P_E2; // edge 2 is the closest
    }

    if (param(0, 0) <= 0 && param(0, 2) >= 1) {
        return PointTriangleDistanceType::P_T0; // vertex 0 is the closest
    } else if (param(0, 1) <= 0 && param(0, 0) >= 1) {
        return PointTriangleDistanceType::P_T1; // vertex 1 is the closest
    } else if (param(0, 2) <= 0 && param(0, 1) >= 1) {
        return PointTriangleDistanceType::P_T2; // vertex 2 is the closest
    } else {
        return PointTriangleDistanceType::P_T;
    }
}

// A more robust implementation of http://geomalgorithms.com/a07-_distance.html
template <typename T>
EdgeEdgeDistanceType edge_edge_distance_type(
    Eigen::ConstRef<Eigen::Vector3<T>> ea0,
    Eigen::ConstRef<Eigen::Vector3<T>> ea1,
    Eigen::ConstRef<Eigen::Vector3<T>> eb0,
    Eigen::ConstRef<Eigen::Vector3<T>> eb1)
{
    // Relative sin² threshold for parallelism: treat edges as parallel when
    // sin²(θ) < PARALLEL_THRESHOLD, scaled by a*c. This avoids misclassifying
    // short nearly-collinear coplanar segments (which an absolute threshold
    // could not handle) and routes such cases to
    // edge_edge_parallel_distance_type.
    //
    // NOTE: u×v cancels for nearly parallel edges, leaving an absolute error of
    // ≈ε‖u‖‖v‖ per component, so sin²(θ) cannot be resolved below ≈ε². The
    // threshold must therefore scale with the precision of T: 2.5e-16 is tuned
    // for double (≈1.13ε), but in single precision it sits ~57× *below* the
    // noise floor, which would make this branch unreachable.
    constexpr T PARALLEL_THRESHOLD = static_cast<T>(
        2.5e-16
        * (static_cast<double>(std::numeric_limits<T>::epsilon())
           / std::numeric_limits<double>::epsilon()));

    const Eigen::Vector3<T> u = ea1 - ea0;
    const Eigen::Vector3<T> v = eb1 - eb0;
    const Eigen::Vector3<T> w = ea0 - eb0;

    const T a = u.squaredNorm(); // always ≥ 0
    const T b = u.dot(v);
    const T c = v.squaredNorm(); // always ≥ 0
    const T d = u.dot(w);
    const T e = v.dot(w);
    const T D = a * c - b * b; // always ≥ 0

    // Degenerate cases should not happen in practice, but we handle them
    if (a == 0 && c == 0) {
        return EdgeEdgeDistanceType::EA0_EB0;
    } else if (a == 0) {
        return EdgeEdgeDistanceType::EA0_EB;
    } else if (c == 0) {
        return EdgeEdgeDistanceType::EA_EB0;
    }

    // Special handling for parallel edges
    const T parallel_tolerance = PARALLEL_THRESHOLD * a * c;
    if (u.cross(v).squaredNorm() < parallel_tolerance) {
        return edge_edge_parallel_distance_type(ea0, ea1, eb0, eb1);
    }

    EdgeEdgeDistanceType default_case = EdgeEdgeDistanceType::EA_EB;

    // compute the line parameters of the two closest points
    const T sN = (b * e - c * d);
    T tN, tD;      // tc = tN / tD
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
        if (tN > 0 && tN < tD
            && u.cross(v).squaredNorm() < parallel_tolerance) {
            // avoid coplanar or nearly parallel EE
            if (sN < D / 2) {
                tN = e;
                tD = c;
                default_case = EdgeEdgeDistanceType::EA0_EB;
            } else {
                tN = e + b;
                tD = c;
                default_case = EdgeEdgeDistanceType::EA1_EB;
            }
        }
        // else default_case stays EdgeEdgeDistanceType::EA_EB
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

template <typename T>
EdgeEdgeDistanceType edge_edge_parallel_distance_type(
    Eigen::ConstRef<Eigen::Vector3<T>> ea0,
    Eigen::ConstRef<Eigen::Vector3<T>> ea1,
    Eigen::ConstRef<Eigen::Vector3<T>> eb0,
    Eigen::ConstRef<Eigen::Vector3<T>> eb1)
{
    const Eigen::Vector3<T> ea = ea1 - ea0;
    const T alpha = (eb0 - ea0).dot(ea) / ea.squaredNorm();
    const T beta = (eb1 - ea0).dot(ea) / ea.squaredNorm();

    uint8_t eac; // 0: EA0, 1: EA1, 2: EA
    uint8_t ebc; // 0: EB0, 1: EB1, 2: EB
    if (alpha < 0) {
        eac = (0 <= beta && beta <= 1) ? 2 : 0;
        ebc = (beta <= alpha) ? 0 : (beta <= 1 ? 1 : 2);
    } else if (alpha > 1) {
        eac = (0 <= beta && beta <= 1) ? 2 : 1;
        ebc = (beta >= alpha) ? 0 : (0 <= beta ? 1 : 2);
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

// clang-format off
template PointEdgeDistanceType point_edge_distance_type<float>(Eigen::ConstRef<VectorMax3f>, Eigen::ConstRef<VectorMax3f>, Eigen::ConstRef<VectorMax3f>);
template PointEdgeDistanceType point_edge_distance_type<double>(Eigen::ConstRef<VectorMax3d>, Eigen::ConstRef<VectorMax3d>, Eigen::ConstRef<VectorMax3d>);
template PointTriangleDistanceType point_triangle_distance_type<float>(Eigen::ConstRef<Eigen::Vector3f>, Eigen::ConstRef<Eigen::Vector3f>, Eigen::ConstRef<Eigen::Vector3f>, Eigen::ConstRef<Eigen::Vector3f>);
template PointTriangleDistanceType point_triangle_distance_type<double>(Eigen::ConstRef<Eigen::Vector3d>, Eigen::ConstRef<Eigen::Vector3d>, Eigen::ConstRef<Eigen::Vector3d>, Eigen::ConstRef<Eigen::Vector3d>);
template EdgeEdgeDistanceType edge_edge_distance_type<float>(Eigen::ConstRef<Eigen::Vector3f>, Eigen::ConstRef<Eigen::Vector3f>, Eigen::ConstRef<Eigen::Vector3f>, Eigen::ConstRef<Eigen::Vector3f>);
template EdgeEdgeDistanceType edge_edge_distance_type<double>(Eigen::ConstRef<Eigen::Vector3d>, Eigen::ConstRef<Eigen::Vector3d>, Eigen::ConstRef<Eigen::Vector3d>, Eigen::ConstRef<Eigen::Vector3d>);
template EdgeEdgeDistanceType edge_edge_parallel_distance_type<float>(Eigen::ConstRef<Eigen::Vector3f>, Eigen::ConstRef<Eigen::Vector3f>, Eigen::ConstRef<Eigen::Vector3f>, Eigen::ConstRef<Eigen::Vector3f>);
template EdgeEdgeDistanceType edge_edge_parallel_distance_type<double>(Eigen::ConstRef<Eigen::Vector3d>, Eigen::ConstRef<Eigen::Vector3d>, Eigen::ConstRef<Eigen::Vector3d>, Eigen::ConstRef<Eigen::Vector3d>);
// clang-format on

} // namespace ipc
