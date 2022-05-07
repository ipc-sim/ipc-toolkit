#pragma once

#include <Eigen/Core>

namespace ipc {

/// @brief Compute the gradient of an edge's length.
/// @param[in] e0 The first vertex of the edge.
/// @param[in] e1 The second vertex of the edge.
/// @param[out] grad The gradient of the edge's length wrt e0, and e1.
template <typename DerivedE0, typename DerivedE1, typename DerivedGrad>
void edge_length_gradient(
    const Eigen::MatrixBase<DerivedE0>& e0,
    const Eigen::MatrixBase<DerivedE1>& e1,
    Eigen::PlainObjectBase<DerivedGrad>& grad)
{
    assert(e0.size() == 2 || e0.size() == 3);
    assert(e1.size() == 2 || e1.size() == 3);
    assert((e1 - e0).norm() != 0);

    // ∇ ‖e₁ - e₀‖
    grad.resize(e0.size() + e1.size());
    grad.head(e0.size()) = (e0 - e1) / (e1 - e0).norm();
    grad.tail(e1.size()) = -grad.head(e0.size());
}

namespace autogen {

    // dA is (9×1) flattened in column-major order
    void triangle_area_gradient(
        double t0_x,
        double t0_y,
        double t0_z,
        double t1_x,
        double t1_y,
        double t1_z,
        double t2_x,
        double t2_y,
        double t2_z,
        double dA[9]);

} // namespace autogen

/// @brief Compute the gradient of the area of a triangle.
/// @param[in] t0 The first vertex of the triangle.
/// @param[in] t1 The second vertex of the triangle.
/// @param[in] t2 The third vertex of the triangle.
/// @param[out] grad The gradient of the triangle's area t0, t1, and t2.
template <
    typename DerivedT0,
    typename DerivedT1,
    typename DerivedT2,
    typename DerivedGrad>
void triangle_area_gradient(
    const Eigen::MatrixBase<DerivedT0>& t0,
    const Eigen::MatrixBase<DerivedT1>& t1,
    const Eigen::MatrixBase<DerivedT2>& t2,
    Eigen::PlainObjectBase<DerivedGrad>& grad)
{
    assert(t0.size() == 3);
    assert(t1.size() == 3);
    assert(t2.size() == 3);

    grad.resize(t0.size() + t1.size() + t2.size());
    autogen::triangle_area_gradient(
        t0[0], t0[1], t0[2], t1[0], t1[1], t1[2], t2[0], t2[1], t2[2],
        grad.data());
}

} // namespace ipc
