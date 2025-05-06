#include "collision.hpp"
#include <ipc/smooth_contact/collisions/smooth_collision.hpp>

namespace ipc {

template <int max_vert>
Collision<max_vert>::Collision(
    const double _weight, const Eigen::SparseVector<double>& _weight_gradient)
    : weight(_weight)
    , weight_gradient(_weight_gradient)
{
}
template <int max_vert>
double Collision<max_vert>::mollifier(
    const Vector<double, -1, Collision<max_vert>::max_size>& positions) const
{
    return 1.0;
}
template <int max_vert>
double Collision<max_vert>::mollifier(
    const Vector<double, -1, Collision<max_vert>::max_size>& positions,
    double eps_x) const
{
    return 1.0;
}
template <int max_vert>
auto Collision<max_vert>::mollifier_gradient(
    const Vector<double, -1, Collision<max_vert>::max_size>& positions) const
    -> Vector<double, -1, max_size>
{
    return Vector<double, -1, Collision<max_vert>::max_size>::Zero(
        positions.size());
}

template <int max_vert>
auto Collision<max_vert>::mollifier_gradient(
    const Vector<double, -1, Collision<max_vert>::max_size>& positions,
    double eps_x) const -> Vector<double, -1, max_size>
{
    return Vector<double, -1, Collision<max_vert>::max_size>::Zero(
        positions.size());
}

template <int max_vert>
auto Collision<max_vert>::mollifier_hessian(
    const Vector<double, -1, Collision<max_vert>::max_size>& positions) const
    -> MatrixMax<double, max_size, max_size>
{
    return MatrixMax<
        double, Collision<max_vert>::max_size, Collision<max_vert>::max_size>::
        Zero(positions.size(), positions.size());
}

template <int max_vert>
auto Collision<max_vert>::mollifier_hessian(
    const Vector<double, -1, Collision<max_vert>::max_size>& positions,
    double eps_x) const -> MatrixMax<double, max_size, max_size>
{
    return MatrixMax<
        double, Collision<max_vert>::max_size, Collision<max_vert>::max_size>::
        Zero(positions.size(), positions.size());
}

template <int max_vert>
Vector12d Collision<max_vert>::mollifier_gradient_wrt_x(
    const Vector<double, -1, Collision<max_vert>::max_size>& rest_positions,
    const Vector<double, -1, Collision<max_vert>::max_size>& positions) const
{
    return Vector12d::Zero(rest_positions.size());
}

template <int max_vert>
Matrix12d Collision<max_vert>::mollifier_gradient_jacobian_wrt_x(
    const Vector<double, -1, Collision<max_vert>::max_size>& rest_positions,
    const Vector<double, -1, Collision<max_vert>::max_size>& positions) const
{
    return Matrix12d::Zero(rest_positions.size(), positions.size());
}

template class Collision<4>;
template class Collision<max_vert_3d>;
template class Collision<max_vert_2d>;

} // namespace ipc
