#include "normal_collision.hpp"

namespace ipc {

NormalCollision::NormalCollision(
    const double _weight, const Eigen::SparseVector<double>& _weight_gradient)
    : weight(_weight)
    , weight_gradient(_weight_gradient)
{
}

double NormalCollision::mollifier(const VectorMax12d& positions) const
{
    return 1.0;
}

double
NormalCollision::mollifier(const VectorMax12d& positions, double eps_x) const
{
    return 1.0;
}

VectorMax12d
NormalCollision::mollifier_gradient(const VectorMax12d& positions) const
{
    return VectorMax12d::Zero(positions.size());
}

VectorMax12d NormalCollision::mollifier_gradient(
    const VectorMax12d& positions, double eps_x) const
{
    return VectorMax12d::Zero(positions.size());
}

MatrixMax12d
NormalCollision::mollifier_hessian(const VectorMax12d& positions) const
{
    return MatrixMax12d::Zero(positions.size(), positions.size());
}

MatrixMax12d NormalCollision::mollifier_hessian(
    const VectorMax12d& positions, double eps_x) const
{
    return MatrixMax12d::Zero(positions.size(), positions.size());
}

Vector12d NormalCollision::mollifier_gradient_wrt_x(
    const VectorMax12d& rest_positions, const VectorMax12d& positions) const
{
    return Vector12d::Zero(rest_positions.size());
}

Matrix12d NormalCollision::mollifier_gradient_jacobian_wrt_x(
    const VectorMax12d& rest_positions, const VectorMax12d& positions) const
{
    return Matrix12d::Zero(rest_positions.size(), positions.size());
}

} // namespace ipc
