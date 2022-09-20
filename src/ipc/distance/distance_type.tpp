#pragma once

#include "distance_type.hpp"

#include <ipc/utils/eigen_ext.hpp>
#include <ipc/utils/logger.hpp>

#include <Eigen/Geometry>

namespace ipc {

template <typename DerivedP, typename DerivedE0, typename DerivedE1>
PointEdgeDistanceType point_edge_distance_type(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedE0>& e0,
    const Eigen::MatrixBase<DerivedE1>& e1)
{
    assert(p.size() == 2 || p.size() == 3);
    assert(e0.size() == 2 || e0.size() == 3);
    assert(e1.size() == 2 || e1.size() == 3);

    const auto e = e1 - e0;
    const auto e_length_sqr = e.squaredNorm();
    if (e_length_sqr == 0) {
        logger().warn("Degenerate edge in point_edge_distance_type!");
        return PointEdgeDistanceType::P_E0; // WARNING: use arbitrary end-point
    }
    auto ratio = e.dot(p - e0) / e_length_sqr;
    if (ratio < 0) {
        return PointEdgeDistanceType::P_E0; // PP (p-e0)
    } else if (ratio > 1) {
        return PointEdgeDistanceType::P_E1; // PP (p-e1)
    } else {
        return PointEdgeDistanceType::P_E; // PE
    }
}

template <
    typename DerivedP,
    typename DerivedT0,
    typename DerivedT1,
    typename DerivedT2>
PointTriangleDistanceType point_triangle_distance_type(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedT0>& t0,
    const Eigen::MatrixBase<DerivedT1>& t1,
    const Eigen::MatrixBase<DerivedT2>& t2)
{
    typedef typename DerivedP::Scalar T;

    assert(p.size() == 3);
    assert(t0.size() == 3);
    assert(t1.size() == 3);
    assert(t2.size() == 3);

    Eigen::Matrix<T, 2, 3> basis;
    basis.row(0) = t1 - t0;
    basis.row(1) = t2 - t0;

    auto normal = cross(basis.row(0), basis.row(1));

    Eigen::Matrix<T, 2, 3> param;

    basis.row(1) = cross(basis.row(0), normal);
    param.col(0) =
        (basis * basis.transpose()).ldlt().solve(basis * Vector3<T>(p - t0));
    if (param(0, 0) > 0.0 && param(0, 0) < 1.0 && param(1, 0) >= 0.0) {
        return PointTriangleDistanceType::P_E0; // edge 0 is the closest
    } else {
        basis.row(0) = t2 - t1;
        basis.row(1) = cross(basis.row(0), normal);
        param.col(1) = (basis * basis.transpose())
                           .ldlt()
                           .solve(basis * Vector3<T>(p - t1));
        if (param(0, 1) > 0.0 && param(0, 1) < 1.0 && param(1, 1) >= 0.0) {
            return PointTriangleDistanceType::P_E1; // edge 1 is the closest
        } else {
            basis.row(0) = t0 - t2;

            basis.row(1) = cross(basis.row(0), normal);
            param.col(2) = (basis * basis.transpose())
                               .ldlt()
                               .solve(basis * Vector3<T>(p - t2));
            if (param(0, 2) > 0.0 && param(0, 2) < 1.0 && param(1, 2) >= 0.0) {
                return PointTriangleDistanceType::P_E2; // edge 2 is the closest
            } else {
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
        }
    }
}

// A more robust implementation of http://geomalgorithms.com/a07-_distance.html
template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1>
EdgeEdgeDistanceType edge_edge_distance_type(
    const Eigen::MatrixBase<DerivedEA0>& ea0,
    const Eigen::MatrixBase<DerivedEA1>& ea1,
    const Eigen::MatrixBase<DerivedEB0>& eb0,
    const Eigen::MatrixBase<DerivedEB1>& eb1)
{
    assert(ea0.size() == 3);
    assert(ea1.size() == 3);
    assert(eb0.size() == 3);
    assert(eb1.size() == 3);

    auto u = ea1 - ea0;
    auto v = eb1 - eb0;
    auto w = ea0 - eb0;

    auto a = u.squaredNorm(); // always >= 0
    auto b = u.dot(v);
    auto c = v.squaredNorm(); // always >= 0
    auto d = u.dot(w);
    auto e = v.dot(w);
    auto D = a * c - b * b; // always >= 0
    auto tD = D;            // tc = tN / tD, default tD = D >= 0

    EdgeEdgeDistanceType defaultCase = EdgeEdgeDistanceType::EA_EB;

    // compute the line parameters of the two closest points
    auto sN = (b * e - c * d);
    decltype(sN) tN;
    if (sN <= 0.0) { // sc < 0 => the s=0 edge is visible
        tN = e;
        tD = c;
        defaultCase = EdgeEdgeDistanceType::EA0_EB;
    } else if (sN >= D) { // sc > 1  => the s=1 edge is visible
        tN = e + b;
        tD = c;
        defaultCase = EdgeEdgeDistanceType::EA1_EB;
    } else {
        tN = (a * e - b * d);
        if (tN > 0.0 && tN < tD
            && (cross(u, v).squaredNorm() < 1.0e-20 * a * c)) {
            // avoid nearly parallel EE
            if (sN < D / 2) {
                tN = e;
                tD = c;
                defaultCase = EdgeEdgeDistanceType::EA0_EB;
            } else {
                tN = e + b;
                tD = c;
                defaultCase = EdgeEdgeDistanceType::EA1_EB;
            }
        }
        // else defaultCase stays EdgeEdgeDistanceType::EA_EB
    }

    if (tN <= 0.0) { // tc < 0 => the t=0 edge is visible
        // recompute sc for this edge
        if (-d <= 0.0) {
            return EdgeEdgeDistanceType::EA0_EB0;
        } else if (-d >= a) {
            return EdgeEdgeDistanceType::EA1_EB0;
        } else {
            return EdgeEdgeDistanceType::EA_EB0;
        }
    } else if (tN >= tD) { // tc > 1  => the t=1 edge is visible
        // recompute sc for this edge
        if ((-d + b) <= 0.0) {
            return EdgeEdgeDistanceType::EA0_EB1;
        } else if ((-d + b) >= a) {
            return EdgeEdgeDistanceType::EA1_EB1;
        } else {
            return EdgeEdgeDistanceType::EA_EB1;
        }
    }

    return defaultCase;
}

} // namespace ipc
