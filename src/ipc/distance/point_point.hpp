#pragma once

#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

/// @brief Compute the distance between two points.
/// @note The distance is actually squared distance.
/// @param p0 The first point.
/// @param p1 The second point.
/// @return The distance between p0 and p1.
template <typename T>
inline T point_point_distance(
    Eigen::ConstRef<VectorMax3<T>> p0, Eigen::ConstRef<VectorMax3<T>> p1)
{
    return (p1 - p0).squaredNorm();
}

/// @brief Compute the gradient of the distance between two points.
/// @note The distance is actually squared distance.
/// @param p0 The first point.
/// @param p1 The second point.
/// @return The computed gradient.
template <typename T>
inline VectorMax6<T> point_point_distance_gradient(
    Eigen::ConstRef<VectorMax3<T>> p0, Eigen::ConstRef<VectorMax3<T>> p1)
{
    const int dim = p0.size();
    assert(p1.size() == dim);

    VectorMax6<T> grad(2 * dim);
    grad.head(dim) = T(2) * (p0 - p1);
    grad.tail(dim) = -grad.head(dim);
    return grad;
}

/// @brief Compute the hessian of the distance between two points.
/// @note The distance is actually squared distance.
/// @param p0 The first point.
/// @param p1 The second point.
/// @return The computed hessian.
template <typename T>
inline MatrixMax6<T> point_point_distance_hessian(
    Eigen::ConstRef<VectorMax3<T>> p0, Eigen::ConstRef<VectorMax3<T>> p1)
{
    const int dim = p0.size();
    assert(p1.size() == dim);

    MatrixMax6<T> hess(2 * dim, 2 * dim);
    hess.setZero();
    hess.diagonal().setConstant(T(2));
    for (int i = 0; i < dim; i++) {
        hess(i, i + dim) = hess(i + dim, i) = T(-2);
    }
    return hess;
}

// --- EigenExpression wrappers ---

template <typename DerivedP0, typename DerivedP1>
inline auto point_point_distance(
    const Eigen::MatrixBase<DerivedP0>& p0,
    const Eigen::MatrixBase<DerivedP1>& p1) -> typename DerivedP0::Scalar
{
    using T = typename DerivedP0::Scalar;
    return point_point_distance(
        Eigen::Ref<const VectorMax3<T>>(p0),
        Eigen::Ref<const VectorMax3<T>>(p1));
}

template <typename DerivedP0, typename DerivedP1>
inline auto point_point_distance_gradient(
    const Eigen::MatrixBase<DerivedP0>& p0,
    const Eigen::MatrixBase<DerivedP1>& p1)
    -> VectorMax6<typename DerivedP0::Scalar>
{
    using T = typename DerivedP0::Scalar;
    return point_point_distance_gradient(
        Eigen::Ref<const VectorMax3<T>>(p0),
        Eigen::Ref<const VectorMax3<T>>(p1));
}

template <typename DerivedP0, typename DerivedP1>
inline auto point_point_distance_hessian(
    const Eigen::MatrixBase<DerivedP0>& p0,
    const Eigen::MatrixBase<DerivedP1>& p1)
    -> MatrixMax6<typename DerivedP0::Scalar>
{
    using T = typename DerivedP0::Scalar;
    return point_point_distance_hessian(
        Eigen::Ref<const VectorMax3<T>>(p0),
        Eigen::Ref<const VectorMax3<T>>(p1));
}

} // namespace ipc
