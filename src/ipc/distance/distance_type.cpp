#include "distance_type.hpp"

#include <ipc/utils/eigen_ext.hpp>
#include <ipc/utils/logger.hpp>

#include <Eigen/Geometry>

namespace ipc {

PointEdgeDistanceType point_edge_distance_type(
    const Eigen::Ref<const VectorMax3d>& p,
    const Eigen::Ref<const VectorMax3d>& e0,
    const Eigen::Ref<const VectorMax3d>& e1)
{
    assert(p.size() == 2 || p.size() == 3);
    assert(e0.size() == 2 || e0.size() == 3);
    assert(e1.size() == 2 || e1.size() == 3);

    const VectorMax3d e = e1 - e0;
    const double e_length_sqr = e.squaredNorm();
    if (e_length_sqr == 0) {
        logger().warn("Degenerate edge in point_edge_distance_type!");
        return PointEdgeDistanceType::P_E0; // WARNING: use arbitrary end-point
    }

    const double ratio = e.dot(p - e0) / e_length_sqr;
    if (ratio < 0) {
        return PointEdgeDistanceType::P_E0; // PP (p-e0)
    } else if (ratio > 1) {
        return PointEdgeDistanceType::P_E1; // PP (p-e1)
    } else {
        return PointEdgeDistanceType::P_E; // PE
    }
}

PointTriangleDistanceType point_triangle_distance_type(
    const Eigen::Ref<const Eigen::Vector3d>& p,
    const Eigen::Ref<const Eigen::Vector3d>& t0,
    const Eigen::Ref<const Eigen::Vector3d>& t1,
    const Eigen::Ref<const Eigen::Vector3d>& t2)
{
    const Eigen::Vector3d normal = (t1 - t0).cross(t2 - t0);

    Eigen::Matrix<double, 2, 3> basis, param;

    basis.row(0) = t1 - t0;
    basis.row(1) = basis.row(0).cross(normal);
    param.col(0) = (basis * basis.transpose()).ldlt().solve(basis * (p - t0));
    if (param(0, 0) > 0.0 && param(0, 0) < 1.0 && param(1, 0) >= 0.0) {
        return PointTriangleDistanceType::P_E0; // edge 0 is the closest
    }

    basis.row(0) = t2 - t1;
    basis.row(1) = basis.row(0).cross(normal);
    param.col(1) = (basis * basis.transpose()).ldlt().solve(basis * (p - t1));
    if (param(0, 1) > 0.0 && param(0, 1) < 1.0 && param(1, 1) >= 0.0) {
        return PointTriangleDistanceType::P_E1; // edge 1 is the closest
    }

    basis.row(0) = t0 - t2;
    basis.row(1) = basis.row(0).cross(normal);
    param.col(2) = (basis * basis.transpose()).ldlt().solve(basis * (p - t2));
    if (param(0, 2) > 0.0 && param(0, 2) < 1.0 && param(1, 2) >= 0.0) {
        return PointTriangleDistanceType::P_E2; // edge 2 is the closest
    }

    if (param(0, 0) <= 0.0 && param(0, 2) >= 1.0) {
        // vertex 0 is the closest
        return PointTriangleDistanceType::P_T0;
    } else if (param(0, 1) <= 0.0 && param(0, 0) >= 1.0) {
        // vertex 1 is the closest
        return PointTriangleDistanceType::P_T1;
    } else if (param(0, 2) <= 0.0 && param(0, 1) >= 1.0) {
        // vertex 2 is the closest
        return PointTriangleDistanceType::P_T2;
    } else {
        return PointTriangleDistanceType::P_T;
    }
}

// A more robust implementation of http://geomalgorithms.com/a07-_distance.html
EdgeEdgeDistanceType edge_edge_distance_type(
    const Eigen::Ref<const Eigen::Vector3d>& ea0,
    const Eigen::Ref<const Eigen::Vector3d>& ea1,
    const Eigen::Ref<const Eigen::Vector3d>& eb0,
    const Eigen::Ref<const Eigen::Vector3d>& eb1)
{
    constexpr double PARALLEL_THRESHOLD = 1.0e-20;

    const Eigen::Vector3d u = ea1 - ea0;
    const Eigen::Vector3d v = eb1 - eb0;
    const Eigen::Vector3d w = ea0 - eb0;

    const double a = u.squaredNorm(); // always ≥ 0
    const double b = u.dot(v);
    const double c = v.squaredNorm(); // always ≥ 0
    const double d = u.dot(w);
    const double e = v.dot(w);
    const double D = a * c - b * b; // always ≥ 0

    // Degenerate cases should not happen in practice, but we handle them
    if (a == 0.0 && c == 0.0) {
        return EdgeEdgeDistanceType::EA0_EB0;
    } else if (a == 0.0) {
        return EdgeEdgeDistanceType::EA0_EB;
    } else if (c == 0.0) {
        return EdgeEdgeDistanceType::EA_EB0;
    }

    // Special handling for parallel edges
    const double parallel_tolerance = PARALLEL_THRESHOLD * std::max(1.0, a * c);
    if (u.cross(v).squaredNorm() < parallel_tolerance) {
        return edge_edge_parallel_distance_type(ea0, ea1, eb0, eb1);
    }

    EdgeEdgeDistanceType default_case = EdgeEdgeDistanceType::EA_EB;

    // compute the line parameters of the two closest points
    const double sN = (b * e - c * d);
    double tN, tD;   // tc = tN / tD
    if (sN <= 0.0) { // sc < 0 ⟹ the s=0 edge is visible
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
        if (tN > 0.0 && tN < tD
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

    if (tN <= 0.0) { // tc < 0 ⟹ the t=0 edge is visible
        // recompute sc for this edge
        if (-d <= 0.0) {
            return EdgeEdgeDistanceType::EA0_EB0;
        } else if (-d >= a) {
            return EdgeEdgeDistanceType::EA1_EB0;
        } else {
            return EdgeEdgeDistanceType::EA_EB0;
        }
    } else if (tN >= tD) { // tc > 1 ⟹ the t=1 edge is visible
        // recompute sc for this edge
        if ((-d + b) <= 0.0) {
            return EdgeEdgeDistanceType::EA0_EB1;
        } else if ((-d + b) >= a) {
            return EdgeEdgeDistanceType::EA1_EB1;
        } else {
            return EdgeEdgeDistanceType::EA_EB1;
        }
    }

    return default_case;
}

EdgeEdgeDistanceType edge_edge_parallel_distance_type(
    const Eigen::Ref<const Eigen::Vector3d>& ea0,
    const Eigen::Ref<const Eigen::Vector3d>& ea1,
    const Eigen::Ref<const Eigen::Vector3d>& eb0,
    const Eigen::Ref<const Eigen::Vector3d>& eb1)
{
    const Eigen::Vector3d ea = ea1 - ea0;
    const double alpha = (eb0 - ea0).dot(ea) / ea.squaredNorm();
    const double beta = (eb1 - ea0).dot(ea) / ea.squaredNorm();

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

} // namespace ipc
