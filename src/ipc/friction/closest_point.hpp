#pragma once

#include <ipc/utils/eigen_ext.hpp>

#include <Eigen/Core>
#include <Eigen/Cholesky>

namespace ipc {

///////////////////////////////////////////////////////////////////////////////
// Point - Edge

template <
    typename DerivedP,
    typename DerivedE0,
    typename DerivedE1,
    typename T = typename DerivedP::Scalar>
inline T point_edge_closest_point(
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
inline Vector2<T> edge_edge_closest_point(
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

    Vector2<T> rhs;
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
inline Vector2<T> point_triangle_closest_point(
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
    basis.row(0) = RowVector3<T>(t1 - t0); // edge 0
    basis.row(1) = RowVector3<T>(t2 - t0); // edge 1
    Matrix2<T> A = basis * basis.transpose();
    Vector2<T> b = basis * Vector3<T>(p - t0);
    Vector2<T> x = A.ldlt().solve(b);
    assert((A * x - b).norm() < 1e-10);
    return x;
}

///////////////////////////////////////////////////////////////////////////////

namespace autogen {
    // J is (6×1) flattened in column-major order
    void point_edge_closest_point_2D_jacobian(
        double p_x,
        double p_y,
        double e0_x,
        double e0_y,
        double e1_x,
        double e1_y,
        double J[6]);

    // J is (9×1) flattened in column-major order
    void point_edge_closest_point_3D_jacobian(
        double p_x,
        double p_y,
        double p_z,
        double e0_x,
        double e0_y,
        double e0_z,
        double e1_x,
        double e1_y,
        double e1_z,
        double J[9]);

    // J is (2×12) flattened in column-major order
    void edge_edge_closest_point_jacobian(
        double ea0_x,
        double ea0_y,
        double ea0_z,
        double ea1_x,
        double ea1_y,
        double ea1_z,
        double eb0_x,
        double eb0_y,
        double eb0_z,
        double eb1_x,
        double eb1_y,
        double eb1_z,
        double J[24]);

    // J is (2×12) flattened in column-major order
    void point_triangle_closest_point_jacobian(
        double p_x,
        double p_y,
        double p_z,
        double t0_x,
        double t0_y,
        double t0_z,
        double t1_x,
        double t1_y,
        double t1_z,
        double t2_x,
        double t2_y,
        double t2_z,
        double J[24]);
} // namespace autogen

template <typename DerivedP, typename DerivedE0, typename DerivedE1>
inline VectorMax9d point_edge_closest_point_jacobian(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedE0>& e0,
    const Eigen::MatrixBase<DerivedE1>& e1)
{
    if (p.size() == 2) {
        assert(e0.size() == 2 && e1.size() == 2);

        Eigen::Matrix<double, 6, 1> J;

        autogen::point_edge_closest_point_2D_jacobian(
            p[0], p[1], e0[0], e0[1], e1[0], e1[1], J.data());

        return J;
    } else {
        assert(p.size() == 3 && e0.size() == 3 && e1.size() == 3);

        Eigen::Matrix<double, 9, 1> J;

        autogen::point_edge_closest_point_3D_jacobian(
            p[0], p[1], p[2], e0[0], e0[1], e0[2], e1[0], e1[1], e1[2],
            J.data());

        return J;
    }
}

template <
    typename DerivedEA0,
    typename DerivedEA1,
    typename DerivedEB0,
    typename DerivedEB1>
inline Eigen::Matrix<double, 2, 12> edge_edge_closest_point_jacobian(
    const Eigen::MatrixBase<DerivedEA0>& ea0,
    const Eigen::MatrixBase<DerivedEA1>& ea1,
    const Eigen::MatrixBase<DerivedEB0>& eb0,
    const Eigen::MatrixBase<DerivedEB1>& eb1)
{
    assert(ea0.size() == 3 && ea1.size() == 3);
    assert(eb0.size() == 3 && eb1.size() == 3);

    Eigen::Matrix<double, 2, 12> J;

    autogen::edge_edge_closest_point_jacobian(
        ea0[0], ea0[1], ea0[2], ea1[0], ea1[1], ea1[2], eb0[0], eb0[1], eb0[2],
        eb1[0], eb1[1], eb1[2], J.data());

    return J;
}

template <
    typename DerivedP,
    typename DerivedT0,
    typename DerivedT1,
    typename DerivedT2>
inline Eigen::Matrix<double, 2, 12> point_triangle_closest_point_jacobian(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedT0>& t0,
    const Eigen::MatrixBase<DerivedT1>& t1,
    const Eigen::MatrixBase<DerivedT2>& t2)
{
    assert(p.size() == 3 && t0.size() == 3 && t1.size() == 3 && t2.size() == 3);

    Eigen::Matrix<double, 2, 12> J;
    autogen::point_triangle_closest_point_jacobian(
        p[0], p[1], p[2], t0[0], t0[1], t0[2], t1[0], t1[1], t1[2], t2[0],
        t2[1], t2[2], J.data());

    return J;
}

} // namespace ipc
