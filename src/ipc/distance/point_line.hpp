#pragma once

#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

// Symbolically generated derivatives
namespace autogen {
    template <typename T>
    void point_line_distance_gradient_2D(
        T v01, T v02, T v11, T v12, T v21, T v22, T g[6]);
    template <typename T>
    void point_line_distance_gradient_3D(
        T v01, T v02, T v03, T v11, T v12, T v13, T v21, T v22, T v23, T g[9]);
    template <typename T>
    void point_line_distance_hessian_2D(
        T v01, T v02, T v11, T v12, T v21, T v22, T H[36]);
    template <typename T>
    void point_line_distance_hessian_3D(
        T v01, T v02, T v03, T v11, T v12, T v13, T v21, T v22, T v23, T H[81]);
} // namespace autogen

/// @brief Compute the distance between a point and line in 2D or 3D.
/// @note The distance is actually squared distance.
/// @param p The point.
/// @param e0 The first vertex of the edge defining the line.
/// @param e1 The second vertex of the edge defining the line.
/// @return The distance between the point and line.
template <typename T>
inline T point_line_distance(
    Eigen::ConstRef<VectorMax3<T>> p,
    Eigen::ConstRef<VectorMax3<T>> e0,
    Eigen::ConstRef<VectorMax3<T>> e1)
{
    assert(p.size() == 2 || p.size() == 3);
    assert(e0.size() == 2 || e0.size() == 3);
    assert(e1.size() == 2 || e1.size() == 3);

    if (p.size() == 2) {
        const Eigen::Vector<T, 2> e = e1 - e0;
        const T numerator =
            (e[1] * p[0] - e[0] * p[1] + e1[0] * e0[1] - e1[1] * e0[0]);
        return numerator * numerator / e.squaredNorm();
    } else {
        const Eigen::Vector3<T> p_to_e0 = e0 - p;
        const Eigen::Vector3<T> p_to_e1 = e1 - p;
        return p_to_e0.cross(p_to_e1).squaredNorm() / (e1 - e0).squaredNorm();
    }
}

/// @brief Compute the gradient of the distance between a point and line.
/// @note The distance is actually squared distance.
/// @param p The point.
/// @param e0 The first vertex of the edge defining the line.
/// @param e1 The second vertex of the edge defining the line.
/// @return The gradient of the distance wrt p, e0, and e1.
template <typename T>
inline VectorMax9<T> point_line_distance_gradient(
    Eigen::ConstRef<VectorMax3<T>> p,
    Eigen::ConstRef<VectorMax3<T>> e0,
    Eigen::ConstRef<VectorMax3<T>> e1)
{
    const int dim = p.size();
    assert(e0.size() == dim);
    assert(e1.size() == dim);

    VectorMax9<T> grad(3 * dim);
    if (dim == 2) {
        autogen::point_line_distance_gradient_2D(
            p[0], p[1], e0[0], e0[1], e1[0], e1[1], grad.data());
    } else {
        autogen::point_line_distance_gradient_3D(
            p[0], p[1], p[2], e0[0], e0[1], e0[2], e1[0], e1[1], e1[2],
            grad.data());
    }
    return grad;
}

/// @brief Compute the hessian of the distance between a point and line.
/// @note The distance is actually squared distance.
/// @param p The point.
/// @param e0 The first vertex of the edge defining the line.
/// @param e1 The second vertex of the edge defining the line.
/// @return The hessian of the distance wrt p, e0, and e1.
template <typename T>
inline MatrixMax9<T> point_line_distance_hessian(
    Eigen::ConstRef<VectorMax3<T>> p,
    Eigen::ConstRef<VectorMax3<T>> e0,
    Eigen::ConstRef<VectorMax3<T>> e1)
{
    const int dim = p.size();
    assert(e0.size() == dim);
    assert(e1.size() == dim);

    MatrixMax9<T> hess(3 * dim, 3 * dim);
    if (dim == 2) {
        autogen::point_line_distance_hessian_2D(
            p[0], p[1], e0[0], e0[1], e1[0], e1[1], hess.data());
    } else {
        autogen::point_line_distance_hessian_3D(
            p[0], p[1], p[2], e0[0], e0[1], e0[2], e1[0], e1[1], e1[2],
            hess.data());
    }
    return hess;
}

// --- EigenExpression wrappers ---

template <
    EigenExpression DerivedP,
    EigenExpression DerivedE0,
    EigenExpression DerivedE1>
inline auto point_line_distance(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedE0>& e0,
    const Eigen::MatrixBase<DerivedE1>& e1) -> typename DerivedP::Scalar
{
    using T = typename DerivedP::Scalar;
    return point_line_distance(
        Eigen::Ref<const VectorMax3<T>>(p), Eigen::Ref<const VectorMax3<T>>(e0),
        Eigen::Ref<const VectorMax3<T>>(e1));
}

template <
    EigenExpression DerivedP,
    EigenExpression DerivedE0,
    EigenExpression DerivedE1>
inline auto point_line_distance_gradient(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedE0>& e0,
    const Eigen::MatrixBase<DerivedE1>& e1)
    -> VectorMax9<typename DerivedP::Scalar>
{
    using T = typename DerivedP::Scalar;
    return point_line_distance_gradient(
        Eigen::Ref<const VectorMax3<T>>(p), Eigen::Ref<const VectorMax3<T>>(e0),
        Eigen::Ref<const VectorMax3<T>>(e1));
}

template <
    EigenExpression DerivedP,
    EigenExpression DerivedE0,
    EigenExpression DerivedE1>
inline auto point_line_distance_hessian(
    const Eigen::MatrixBase<DerivedP>& p,
    const Eigen::MatrixBase<DerivedE0>& e0,
    const Eigen::MatrixBase<DerivedE1>& e1)
    -> MatrixMax9<typename DerivedP::Scalar>
{
    using T = typename DerivedP::Scalar;
    return point_line_distance_hessian(
        Eigen::Ref<const VectorMax3<T>>(p), Eigen::Ref<const VectorMax3<T>>(e0),
        Eigen::Ref<const VectorMax3<T>>(e1));
}

} // namespace ipc
