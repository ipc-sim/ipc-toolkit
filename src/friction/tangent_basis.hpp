#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <utils/eigen_ext.hpp>

namespace ipc {

// Point - Triangle

/// Compute a basis for the space tangent to the point-triangle pair.
template <
    typename DerivedP,
    typename DerivedT0,
    typename DerivedT1,
    typename DerivedT2>
inline auto point_triangle_tangent_basis(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedT0>& t0,
    const Eigen::MatrixBase<DerivedT1>& t1,
    const Eigen::MatrixBase<DerivedT2>& t2)
{
    assert(p.size() == 3);
    assert(t0.size() == 3);
    assert(t1.size() == 3);
    assert(t2.size() == 3);

    typedef typename DerivedP::Scalar T;

    Eigen::Matrix<T, 3, 2> basis;

    auto e0 = t1 - t0;
    // The first basis vector is along first edge of the triangle.
    basis.col(0) = e0.normalized();
    // The second basis vector is orthogonal to the first and the triangle
    // normal.
    auto normal = Eigen::cross(e0, t2 - t0);
    assert(normal.norm() != 0);
    basis.col(1) = Eigen::cross(normal, e0).normalized();

    return basis;
}

// Edge - Edge

/// Compute a basis for the space tangent to the edge-edge pair.
template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1>
inline auto edge_edge_tangent_basis(
    const Eigen::MatrixBase<DerivedEA0>& ea0,
    const Eigen::MatrixBase<DerivedEA1>& ea1,
    const Eigen::MatrixBase<DerivedEB0>& eb0,
    const Eigen::MatrixBase<DerivedEB1>& eb1)
{
    assert(ea0.size() == 3);
    assert(ea1.size() == 3);
    assert(eb0.size() == 3);
    assert(eb1.size() == 3);

    typedef typename DerivedEA0::Scalar T;

    Eigen::Matrix<T, 3, 2> basis;

    auto ea = ea1 - ea0; // Edge A direction
    // The first basis vector is along edge A.
    basis.col(0) = ea0.normalized();
    // The second basis vector is orthogonal to the first and the edge-edge
    // normal.
    auto normal = Eigen::cross(ea, eb1 - eb0);
    // The normal will be zero if the edges are parallel (i.e. coplanar).
    assert(normal.norm() != 0);
    basis.col(1) = Eigen::cross(normal, ea).normalized();

    return basis;
}

// Point - Edge

template <typename DerivedP, typename DerivedE0, typename DerivedE1>
inline auto point_edge_tangent_basis(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedE0>& e0,
    const Eigen::MatrixBase<DerivedE1>& e1)
{
    assert(p.size() == 3);
    assert(e0.size() == 3);
    assert(e1.size() == 3);

    typedef typename DerivedP::Scalar T;

    Eigen::Matrix<T, 3, 2> basis;

    auto e = e1 - e0;
    basis.col(0) = e.normalized();
    basis.col(1) = Eigen::cross(e, p - e0).normalized();

    return basis;
}

// Point - Point

template <typename DerivedP0, typename DerivedP1>
inline auto point_point_tangent_basis(
    const Eigen::MatrixBase<DerivedP0>& p0,
    const Eigen::MatrixBase<DerivedP1>& p1)
{
    assert(p0.size() == 3);
    assert(p1.size() == 3);

    typedef typename DerivedP0::Scalar T;

    Eigen::Matrix<T, 3, 2> basis;

    auto p0_to_p1 = (p1 - p0).transpose();

    auto cross_x = Eigen::cross(Eigen::Matrix<T, 1, 3>::UnitX(), p0_to_p1);
    auto cross_y = Eigen::cross(Eigen::Matrix<T, 1, 3>::UnitY(), p0_to_p1);

    if (cross_x.squaredNorm() > cross_y.squaredNorm()) {
        basis.col(0) = cross_x.normalized().transpose();
        basis.col(1) = Eigen::cross(p0_to_p1, cross_x).normalized().transpose();
    } else {
        basis.col(0) = cross_y.normalized().transpose();
        basis.col(1) = Eigen::cross(p0_to_p1, cross_y).normalized().transpose();
    }

    return basis;
}

} // namespace ipc
