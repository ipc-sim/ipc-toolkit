#include "smooth_contact_potential.hpp"

namespace ipc {

// -- Single collision methods -------------------------------------------------

template <class TCollisions>
double SmoothContactPotential<TCollisions>::operator()(
    const Collision& collision, const VectorMax12d& positions) const
{
    return collision.weight * collision(positions, params);
}

template <class TCollisions>
VectorMax12d SmoothContactPotential<TCollisions>::gradient(
    const Collision& collision, const VectorMax12d& positions) const
{
    return collision.weight * collision.gradient(positions, params);
}

template <class TCollisions>
MatrixMax12d SmoothContactPotential<TCollisions>::hessian(
    const Collision& collision,
    const VectorMax12d& positions,
    const bool project_hessian_to_psd) const
{
    return collision.weight * collision.hessian(positions, params, project_hessian_to_psd);
}

template class SmoothContactPotential<VirtualCollisions>;
template class SmoothContactPotential<SmoothCollisions>;

} // namespace ipc