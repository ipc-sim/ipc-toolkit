#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

// Point - Point

template <
    typename DerivedP0,
    typename DerivedP1,
    typename T = typename DerivedP0::Scalar>
inline Eigen::MatrixXX<T, 3, 2> point_point_tangent_basis(
    const Eigen::MatrixBase<DerivedP0>& p0,
    const Eigen::MatrixBase<DerivedP1>& p1)
{
    if (p0.size() == 2) {
        assert(p1.size() == 2);

        Eigen::MatrixXX<T, 3, 2> basis(2, 1);

        auto p0_to_p1 = (p1 - p0).normalized();

        basis(0) = -p0_to_p1(1);
        basis(1) = p0_to_p1(0);

        return basis;
    } else {
        assert(p0.size() == 3 && p1.size() == 3);

        Eigen::MatrixXX<T, 3, 2> basis(3, 2);

        auto p0_to_p1 = p1 - p0;

        Eigen::Vector3<T> cross_x =
            Eigen::cross(Eigen::Vector3<T>::UnitX(), p0_to_p1);
        Eigen::Vector3<T> cross_y =
            Eigen::cross(Eigen::Vector3<T>::UnitY(), p0_to_p1);

        if (cross_x.squaredNorm() > cross_y.squaredNorm()) {
            basis.col(0) = cross_x.normalized();
            basis.col(1) = Eigen::cross(p0_to_p1, cross_x).normalized();
        } else {
            basis.col(0) = cross_y.normalized();
            basis.col(1) = Eigen::cross(p0_to_p1, cross_y).normalized();
        }

        return basis;
    }
}

///////////////////////////////////////////////////////////////////////////////
// Point - Edge

template <
    typename DerivedP,
    typename DerivedE0,
    typename DerivedE1,
    typename T = typename DerivedP::Scalar>
inline Eigen::MatrixXX<T, 3, 2> point_edge_tangent_basis(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedE0>& e0,
    const Eigen::MatrixBase<DerivedE1>& e1)
{
    if (p.size() == 2) {
        assert(e0.size() == 2 && e1.size() == 2);

        Eigen::MatrixXX<T, 3, 2> basis(2, 1);

        basis.col(0) = (e1 - e0).normalized();

        return basis;
    } else {
        assert(p.size() == 3 && e0.size() == 3 && e1.size() == 3);

        Eigen::MatrixXX<T, 3, 2> basis(3, 2);

        auto e = e1 - e0;
        basis.col(0) = e.normalized();
        basis.col(1) = Eigen::cross(e, p - e0).normalized();

        return basis;
    }
}

///////////////////////////////////////////////////////////////////////////////
// Edge - Edge

/// Compute a basis for the space tangent to the edge-edge pair.
template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1,
    typename T = typename DerivedEA0::Scalar>
inline Eigen::Matrix<T, 3, 2> edge_edge_tangent_basis(
    const Eigen::MatrixBase<DerivedEA0>& ea0,
    const Eigen::MatrixBase<DerivedEA1>& ea1,
    const Eigen::MatrixBase<DerivedEB0>& eb0,
    const Eigen::MatrixBase<DerivedEB1>& eb1)
{
    assert(ea0.size() == 3 && ea1.size() == 3);
    assert(eb0.size() == 3 && eb1.size() == 3);

    Eigen::Matrix<T, 3, 2> basis;

    auto ea = ea1 - ea0; // Edge A direction
    // The first basis vector is along edge A.
    basis.col(0) = ea.normalized();
    // The second basis vector is orthogonal to the first and the edge-edge
    // normal.
    auto normal = Eigen::cross(ea, eb1 - eb0);
    // The normal will be zero if the edges are parallel (i.e. coplanar).
    assert(normal.norm() != 0);
    basis.col(1) = Eigen::cross(normal, ea).normalized();

    return basis;
}

///////////////////////////////////////////////////////////////////////////////
// Point - Triangle

/// Compute a basis for the space tangent to the point-triangle pair.
template <
    typename DerivedP,
    typename DerivedT0,
    typename DerivedT1,
    typename DerivedT2,
    typename T = typename DerivedP::Scalar>
inline Eigen::Matrix<T, 3, 2> point_triangle_tangent_basis(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedT0>& t0,
    const Eigen::MatrixBase<DerivedT1>& t1,
    const Eigen::MatrixBase<DerivedT2>& t2)
{
    assert(p.size() == 3 && t0.size() == 3 && t1.size() == 3 && t2.size() == 3);

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

} // namespace ipc
