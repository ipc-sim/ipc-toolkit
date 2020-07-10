#pragma once

#include <Eigen/Core>

namespace ipc {

template <typename DerivedP0, typename DerivedP1>
inline auto point_point_distance(
    const Eigen::MatrixBase<DerivedP0>& p0,
    const Eigen::MatrixBase<DerivedP1>& p1)
{
#ifdef USE_DISTANCE_SQUARED
    return (p1 - p0).squaredNorm();
#else
    return (p1 - p0).norm();
#endif
}

template <typename DerivedP0, typename DerivedP1, typename DerivedGrad>
inline void point_point_distance_gradient(
    const Eigen::MatrixBase<DerivedP0>& p0,
    const Eigen::MatrixBase<DerivedP1>& p1,
    Eigen::MatrixBase<DerivedGrad>& grad)
{
#ifdef USE_DISTANCE_SQUARED
    grad.head(p0.size()) = 2.0 * (p0 - p1);
    grad.tail(p1.size()) = -grad.head(p0.size());
#else
    static_assert(false, "missing non-squared distance definition");
#endif
}

template <typename DerivedP0, typename DerivedP1, typename DerivedHess>
inline void point_point_distance_hessian(
    const Eigen::MatrixBase<DerivedP0>& p0,
    const Eigen::MatrixBase<DerivedP1>& p1,
    Eigen::MatrixBase<DerivedHess>& hess)
{
#ifdef USE_DISTANCE_SQUARED
    int dim = p0.size();
    hess.setZero();
    hess.diagonal().setConstant(2.0);
    for (int i = 0; i < dim; i++) {
        hess(i, i + dim) = hess(i + dim, i) = -2;
    }
#else
    static_assert(false, "missing non-squared distance definition");
#endif
}

} // namespace ipc
