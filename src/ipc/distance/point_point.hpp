#pragma once

#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

/// @brief Compute the distance between two points.
/// @note The distance is actually squared distance.
/// @param[in] p0 The first point.
/// @param[in] p1 The second point.
/// @return The distance between p0 and p1.
template <typename DerivedP0, typename DerivedP1>
inline auto point_point_distance(
    const Eigen::MatrixBase<DerivedP0>& p0,
    const Eigen::MatrixBase<DerivedP1>& p1)
{
    return (p1 - p0).squaredNorm();
}

/// @brief Compute the gradient of the distance between two points.
/// @note The distance is actually squared distance.
/// @param[in] p0 The first point.
/// @param[in] p1 The second point.
/// @param[out] grad The computed gradient.
template <typename DerivedP0, typename DerivedP1, typename DerivedGrad>
inline void point_point_distance_gradient(
    const Eigen::MatrixBase<DerivedP0>& p0,
    const Eigen::MatrixBase<DerivedP1>& p1,
    Eigen::PlainObjectBase<DerivedGrad>& grad)
{
    assert(p0.size() == p1.size());
    grad.resize(p0.size() + p1.size());
    grad.head(p0.size()) = 2.0 * (p0 - p1);
    grad.tail(p1.size()) = -grad.head(p0.size());
}

/// @brief Compute the hessian of the distance between two points.
/// @note The distance is actually squared distance.
/// @param[in] p0 The first point.
/// @param[in] p1 The second point.
/// @param[out] hess The computed hessian.
template <typename DerivedP0, typename DerivedP1, typename DerivedHess>
inline void point_point_distance_hessian(
    const Eigen::MatrixBase<DerivedP0>& p0,
    const Eigen::MatrixBase<DerivedP1>& p1,
    Eigen::PlainObjectBase<DerivedHess>& hess)
{
    int dim = p0.size();
    assert(p1.size() == dim);

    hess.resize(2 * dim, 2 * dim);

    hess.setZero();
    hess.diagonal().setConstant(2.0);
    for (int i = 0; i < dim; i++) {
        hess(i, i + dim) = hess(i + dim, i) = -2;
    }
}

} // namespace ipc
