#pragma once

#include <iostream>

#include <Eigen/Geometry>
#include <igl/barycentric_coordinates.h>

#include <utils/eigen_ext.hpp>

namespace ipc {

template <typename DerivedP, typename DerivedE0, typename DerivedE1>
PointEdgeDistanceType point_edge_distance_type(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedE0>& e0,
    const Eigen::MatrixBase<DerivedE1>& e1)
{
    const auto e = e1 - e0;
    const auto e_length_sqr = e.squaredNorm();
    if (e_length_sqr == 0) {
        return PointEdgeDistanceType::P_E; // PE
    }
    auto ratio = e.dot(p - e0) / e_length_sqr;
    if (ratio <= 0) {
        return PointEdgeDistanceType::P_E0; // PP (p-e0)
    } else if (ratio >= 1) {
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

    // Compute the barycentric coordinates of the projected_point
    const auto normal = Eigen::cross(t1 - t0, t2 - t0);
    const auto projected_p =
        p - (p - t0).dot(normal) / normal.squaredNorm() * normal;

    typedef Eigen::Matrix<T, 1, 3> RowVector3T;

    Eigen::MatrixX<T> coords;
    igl::barycentric_coordinates(
        RowVector3T(projected_p),                          //
        RowVector3T(t0), RowVector3T(t1), RowVector3T(t2), //
        coords);
    T u = coords(0, 0), v = coords(0, 1), w = coords(0, 2);
    assert(abs(u + v + w - 1) < 1e-12);

    // Find the closest point using the barycentric coordinates
    // https://math.stackexchange.com/a/589362

    // Is closest point in the plane inside the trianlge?
    if (u >= 0 && v >= 0 && w >= 0) {
        return PointTriangleDistanceType::P_T;
    }

    // Check if a vertex is the closest point on the triangle
    if (u >= 0 && v < 0 && w < 0) {
        // vertex 0 is the closest
        return PointTriangleDistanceType::P_T0;
    }
    if (u < 0 && v >= 0 && w < 0) {
        // vertex 1 is the closest
        return PointTriangleDistanceType::P_T1;
    }
    if (u < 0 && v < 0 && w >= 0) {
        // vertex 2 is the closest
        return PointTriangleDistanceType::P_T2;
    }

    // Check if an edge is the closest point on the triangle
    if (u >= 0 && v >= 0 && w < 0) {
        // edge 0 is the closest
        PointEdgeDistanceType pe_dtype = point_edge_distance_type(p, t0, t1);
        switch (pe_dtype) {
        case PointEdgeDistanceType::P_E0:
            return PointTriangleDistanceType::P_T0;
        case PointEdgeDistanceType::P_E1:
            return PointTriangleDistanceType::P_T1;
        case PointEdgeDistanceType::P_E:
            return PointTriangleDistanceType::P_E0;
        }
    }
    if (u < 0 && v >= 0 && w >= 0) {
        // edge 1 is the closest
        PointEdgeDistanceType pe_dtype = point_edge_distance_type(p, t1, t2);
        switch (pe_dtype) {
        case PointEdgeDistanceType::P_E0:
            return PointTriangleDistanceType::P_T1;
        case PointEdgeDistanceType::P_E1:
            return PointTriangleDistanceType::P_T2;
        case PointEdgeDistanceType::P_E:
            return PointTriangleDistanceType::P_E1;
        }
    }
    if (u >= 0 && v < 0 && w >= 0) {
        // edge 2 is the closest
        PointEdgeDistanceType pe_dtype = point_edge_distance_type(p, t2, t0);
        switch (pe_dtype) {
        case PointEdgeDistanceType::P_E0:
            return PointTriangleDistanceType::P_T2;
        case PointEdgeDistanceType::P_E1:
            return PointTriangleDistanceType::P_T0;
        case PointEdgeDistanceType::P_E:
            return PointTriangleDistanceType::P_E2;
        }
    }

    // This should never happen because u + v + w = 1.
    throw "point_triangle_distance is not implemented correctly!";
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
            && (Eigen::cross(u, v).squaredNorm() < 1.0e-20 * a * c)) {
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
        // else defaultCase stays as 8
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
