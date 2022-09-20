#pragma once

#include <ipc/utils/eigen_ext.hpp>

#include <Eigen/Geometry>

namespace ipc {

// Point - Point

template <
    typename DerivedP0,
    typename DerivedP1,
    typename T = typename DerivedP0::Scalar>
inline MatrixMax<T, 3, 2> point_point_tangent_basis(
    const Eigen::MatrixBase<DerivedP0>& p0,
    const Eigen::MatrixBase<DerivedP1>& p1)
{
    if (p0.size() == 2) {
        assert(p1.size() == 2);

        MatrixMax<T, 3, 2> basis(2, 1);

        auto p0_to_p1 = (p1 - p0).normalized();

        basis(0) = -p0_to_p1(1);
        basis(1) = p0_to_p1(0);

        return basis;
    } else {
        assert(p0.size() == 3 && p1.size() == 3);

        MatrixMax<T, 3, 2> basis(3, 2);

        auto p0_to_p1 = p1 - p0;

        Vector3<T> cross_x = cross(Vector3<T>::UnitX(), p0_to_p1);
        Vector3<T> cross_y = cross(Vector3<T>::UnitY(), p0_to_p1);

        if (cross_x.squaredNorm() > cross_y.squaredNorm()) {
            basis.col(0) = cross_x.normalized();
            basis.col(1) = cross(p0_to_p1, cross_x).normalized();
        } else {
            basis.col(0) = cross_y.normalized();
            basis.col(1) = cross(p0_to_p1, cross_y).normalized();
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
inline MatrixMax<T, 3, 2> point_edge_tangent_basis(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedE0>& e0,
    const Eigen::MatrixBase<DerivedE1>& e1)
{
    if (p.size() == 2) {
        assert(e0.size() == 2 && e1.size() == 2);

        MatrixMax<T, 3, 2> basis(2, 1);

        basis.col(0) = (e1 - e0).normalized();

        return basis;
    } else {
        assert(p.size() == 3 && e0.size() == 3 && e1.size() == 3);

        MatrixMax<T, 3, 2> basis(3, 2);

        auto e = e1 - e0;
        basis.col(0) = e.normalized();
        basis.col(1) = cross(e, p - e0).normalized();

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
    auto normal = cross(ea, eb1 - eb0);
    // The normal will be zero if the edges are parallel (i.e. coplanar).
    assert(normal.norm() != 0);
    basis.col(1) = cross(normal, ea).normalized();

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
    auto normal = cross(e0, t2 - t0);
    assert(normal.norm() != 0);
    basis.col(1) = cross(normal, e0).normalized();

    return basis;
}

////////////////////////////////////////////////////////////////////////////////
// Gradient

namespace autogen {
    // J is (8×1) flattened in column-major order
    void point_point_tangent_basis_2D_jacobian(
        double p0_x, double p0_y, double p1_x, double p1_y, double J[8]);

    // J is (18×2) flattened in column-major order
    void point_point_tangent_basis_3D_jacobian(
        double p0_x,
        double p0_y,
        double p0_z,
        double p1_x,
        double p1_y,
        double p1_z,
        double J[36]);

    // J is (12×1) flattened in column-major order
    void point_edge_tangent_basis_2D_jacobian(
        double p_x,
        double p_y,
        double e0_x,
        double e0_y,
        double e1_x,
        double e1_y,
        double J[12]);

    // J is (27×2) flattened in column-major order
    void point_edge_tangent_basis_3D_jacobian(
        double p_x,
        double p_y,
        double p_z,
        double e0_x,
        double e0_y,
        double e0_z,
        double e1_x,
        double e1_y,
        double e1_z,
        double J[54]);

    // J is (36×2) flattened in column-major order
    void edge_edge_tangent_basis_jacobian(
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
        double J[72]);

    // J is (36×2) flattened in column-major order
    void point_triangle_tangent_basis_jacobian(
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
        double J[72]);
} // namespace autogen

template <typename DerivedP0, typename DerivedP1>
inline MatrixMax<double, 18, 2> point_point_tangent_basis_jacobian(
    const Eigen::MatrixBase<DerivedP0>& p0,
    const Eigen::MatrixBase<DerivedP1>& p1)
{
    if (p0.size() == 2) {
        assert(p1.size() == 2);

        Eigen::Matrix<double, 8, 1> J;

        autogen::point_point_tangent_basis_2D_jacobian(
            p0[0], p0[1], p1[0], p1[1], J.data());

        return J;
    } else {
        assert(p0.size() == 3 && p1.size() == 3);

        Eigen::Matrix<double, 18, 2> J;

        autogen::point_point_tangent_basis_3D_jacobian(
            p0[0], p0[1], p0[2], p1[0], p1[1], p1[2], J.data());

        return J;
    }
}

template <typename DerivedP, typename DerivedE0, typename DerivedE1>
inline MatrixMax<double, 27, 2> point_edge_tangent_basis_jacobian(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedE0>& e0,
    const Eigen::MatrixBase<DerivedE1>& e1)
{
    if (p.size() == 2) {
        assert(e0.size() == 2 && e1.size() == 2);

        Eigen::Matrix<double, 12, 1> J;

        autogen::point_edge_tangent_basis_2D_jacobian(
            p[0], p[1], e0[0], e0[1], e1[0], e1[1], J.data());

        return J;
    } else {
        assert(p.size() == 3 && e0.size() == 3 && e1.size() == 3);

        Eigen::Matrix<double, 27, 2> J;

        autogen::point_edge_tangent_basis_3D_jacobian(
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
inline Eigen::Matrix<double, 36, 2> edge_edge_tangent_basis_jacobian(
    const Eigen::MatrixBase<DerivedEA0>& ea0,
    const Eigen::MatrixBase<DerivedEA1>& ea1,
    const Eigen::MatrixBase<DerivedEB0>& eb0,
    const Eigen::MatrixBase<DerivedEB1>& eb1)
{
    assert(ea0.size() == 3 && ea1.size() == 3);
    assert(eb0.size() == 3 && eb1.size() == 3);

    Eigen::Matrix<double, 36, 2> J;

    autogen::edge_edge_tangent_basis_jacobian(
        ea0[0], ea0[1], ea0[2], ea1[0], ea1[1], ea1[2], eb0[0], eb0[1], eb0[2],
        eb1[0], eb1[1], eb1[2], J.data());

    return J;
}

template <
    typename DerivedP,
    typename DerivedT0,
    typename DerivedT1,
    typename DerivedT2>
inline Eigen::Matrix<double, 36, 2> point_triangle_tangent_basis_jacobian(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedT0>& t0,
    const Eigen::MatrixBase<DerivedT1>& t1,
    const Eigen::MatrixBase<DerivedT2>& t2)
{
    assert(p.size() == 3 && t0.size() == 3 && t1.size() == 3 && t2.size() == 3);

    Eigen::Matrix<double, 36, 2> J;

    autogen::point_triangle_tangent_basis_jacobian(
        p[0], p[1], p[2], t0[0], t0[1], t0[2], t1[0], t1[1], t1[2], t2[0],
        t2[1], t2[2], J.data());

    return J;
}

} // namespace ipc
