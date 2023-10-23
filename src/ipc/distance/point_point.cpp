#include "point_point.hpp"

namespace ipc {

double point_point_distance(
    const Eigen::Ref<const VectorMax3d>& p0,
    const Eigen::Ref<const VectorMax3d>& p1)
{
    return (p1 - p0).squaredNorm();
}

VectorMax6d point_point_distance_gradient(
    const Eigen::Ref<const VectorMax3d>& p0,
    const Eigen::Ref<const VectorMax3d>& p1)
{
    int dim = p0.size();
    assert(p1.size() == dim);

    VectorMax6d grad(2 * dim);

    grad.head(dim) = 2.0 * (p0 - p1);
    grad.tail(dim) = -grad.head(dim);

    return grad;
}

MatrixMax6d point_point_distance_hessian(
    const Eigen::Ref<const VectorMax3d>& p0,
    const Eigen::Ref<const VectorMax3d>& p1)
{
    int dim = p0.size();
    assert(p1.size() == dim);

    MatrixMax6d hess(2 * dim, 2 * dim);

    hess.setZero();
    hess.diagonal().setConstant(2.0);
    for (int i = 0; i < dim; i++) {
        hess(i, i + dim) = hess(i + dim, i) = -2;
    }

    return hess;
}

} // namespace ipc
