#include "collision.hpp"

namespace ipc {

template <int max_vert>
Collision<max_vert>::Collision(
    const double _weight, const Eigen::SparseVector<double>& _weight_gradient)
    : weight(_weight)
    , weight_gradient(_weight_gradient)
{
}
template <int max_vert>
double Collision<max_vert>::mollifier(const Vector<double, -1, 3*max_vert>& positions) const { return 1.0; }
template <int max_vert>
double Collision<max_vert>::mollifier(const Vector<double, -1, 3*max_vert>& positions, double eps_x) const
{
    return 1.0;
}
template <int max_vert>
Vector<double, -1, 3*max_vert> Collision<max_vert>::mollifier_gradient(const Vector<double, -1, 3*max_vert>& positions) const
{
    return Vector<double, -1, 3*max_vert>::Zero(positions.size());
}

template <int max_vert> Vector<double, -1, 3*max_vert>
Collision<max_vert>::mollifier_gradient(const Vector<double, -1, 3*max_vert>& positions, double eps_x) const
{
    return Vector<double, -1, 3*max_vert>::Zero(positions.size());
}

template <int max_vert>
MatrixMax<double, 3*max_vert, 3*max_vert> Collision<max_vert>::mollifier_hessian(const Vector<double, -1, 3*max_vert>& positions) const
{
    return MatrixMax<double, 3*max_vert, 3*max_vert>::Zero(positions.size(), positions.size());
}

template <int max_vert> MatrixMax<double, 3*max_vert, 3*max_vert>
Collision<max_vert>::mollifier_hessian(const Vector<double, -1, 3*max_vert>& positions, double eps_x) const
{
    return MatrixMax<double, 3*max_vert, 3*max_vert>::Zero(positions.size(), positions.size());
}

template <int max_vert>
Vector12d Collision<max_vert>::mollifier_gradient_wrt_x(
    const Vector<double, -1, 3*max_vert>& rest_positions, const Vector<double, -1, 3*max_vert>& positions) const
{
    return Vector12d::Zero(rest_positions.size());
}

template <int max_vert>
Matrix12d Collision<max_vert>::mollifier_gradient_jacobian_wrt_x(
    const Vector<double, -1, 3*max_vert>& rest_positions, const Vector<double, -1, 3*max_vert>& positions) const
{
    return Matrix12d::Zero(rest_positions.size(), positions.size());
}

template class Collision<4>;
template class Collision<6>;

} // namespace ipc
