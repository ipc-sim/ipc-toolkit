#pragma once

#include <Eigen/Cholesky>
#include <Eigen/Core>

#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

///////////////////////////////////////////////////////////////////////////////
// Point - Edge

template <typename DerivedP, typename DerivedE0, typename DerivedE1>
inline auto point_edge_closest_point(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedE0>& e0,
    const Eigen::MatrixBase<DerivedE1>& e1)
{
    auto e = e1 - e0;
    return (p - e0).dot(e) / e.squaredNorm();
}

///////////////////////////////////////////////////////////////////////////////
// Edge - Edge

/// Compute the barycentric coordinates of the closest points
template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1,
    typename T = typename DerivedEA0::Scalar>
inline Eigen::Vector2<T> edge_edge_closest_point(
    const Eigen::MatrixBase<DerivedEA0>& ea0,
    const Eigen::MatrixBase<DerivedEA1>& ea1,
    const Eigen::MatrixBase<DerivedEB0>& eb0,
    const Eigen::MatrixBase<DerivedEB1>& eb1)
{
    assert(ea0.size() == 3);
    assert(ea1.size() == 3);
    assert(eb0.size() == 3);
    assert(eb1.size() == 3);

    auto eb_to_ea = ea0 - eb0;
    auto ea = ea1 - ea0;
    auto eb = eb1 - eb0;

    Eigen::Matrix<T, 2, 2> coefMtr;
    coefMtr(0, 0) = ea.squaredNorm();
    coefMtr(0, 1) = coefMtr(1, 0) = -eb.dot(ea);
    coefMtr(1, 1) = eb.squaredNorm();

    Eigen::Vector2<T> rhs;
    rhs[0] = -eb_to_ea.dot(ea);
    rhs[1] = eb_to_ea.dot(eb);

    return coefMtr.ldlt().solve(rhs);
}

///////////////////////////////////////////////////////////////////////////////
// Point - Triangle

/// Compute the barycentric coordinates of the closest point on the triangle.
template <
    typename DerivedP,
    typename DerivedT0,
    typename DerivedT1,
    typename DerivedT2,
    typename T = typename DerivedP::Scalar>
inline Eigen::Vector2<T> point_triangle_closest_point(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedT0>& t0,
    const Eigen::MatrixBase<DerivedT1>& t1,
    const Eigen::MatrixBase<DerivedT2>& t2)
{
    assert(p.size() == 3);
    assert(t0.size() == 3);
    assert(t1.size() == 3);
    assert(t2.size() == 3);

    Eigen::Matrix<T, 2, 3> basis;
    basis.row(0) = Eigen::RowVector3<T>(t1 - t0); // edge 0
    basis.row(1) = Eigen::RowVector3<T>(t2 - t0); // edge 1
    Eigen::Matrix2<T> A = basis * basis.transpose();
    Eigen::Vector2<T> b = basis * Eigen::Vector3<T>(p - t0);
    Eigen::Vector2<T> x = A.ldlt().solve(b);
    assert((A * x - b).norm() < 1e-10);
    return x;
}

} // namespace ipc
