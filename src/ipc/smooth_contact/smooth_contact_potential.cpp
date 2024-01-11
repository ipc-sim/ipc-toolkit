#include "smooth_contact_potential.hpp"

namespace ipc {

// -- Single collision methods -------------------------------------------------

template <class TCollisions>
double SmoothContactPotential<TCollisions>::operator()(
    const typename TCollisions::value_type& collision, 
    const Vector<double, -1, SmoothContactPotential<TCollisions>::element_size>& positions) const
{
    return collision.weight * collision(positions, params);
}

template <class TCollisions>
Vector<double, -1, SmoothContactPotential<TCollisions>::element_size> SmoothContactPotential<TCollisions>::gradient(
    const typename TCollisions::value_type& collision, 
    const Vector<double, -1, SmoothContactPotential<TCollisions>::element_size>& positions) const
{
    return collision.weight * collision.gradient(positions, params);
}

template <class TCollisions>
MatrixMax<double, SmoothContactPotential<TCollisions>::element_size, SmoothContactPotential<TCollisions>::element_size> SmoothContactPotential<TCollisions>::hessian(
    const typename TCollisions::value_type& collision,
    const Vector<double, -1, SmoothContactPotential<TCollisions>::element_size>& positions,
    const bool project_hessian_to_psd) const
{
    return collision.weight * collision.hessian(positions, params, project_hessian_to_psd);
}

// template class SmoothContactPotential<VirtualCollisions<4>>;
template class SmoothContactPotential<SmoothCollisions<2, SmoothEdgeEdgeCollision<2>>>;
template class SmoothContactPotential<SmoothCollisions<3, SmoothFaceFaceCollision>>;

} // namespace ipc