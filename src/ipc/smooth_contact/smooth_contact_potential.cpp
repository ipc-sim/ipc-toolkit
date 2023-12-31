#include "smooth_contact_potential.hpp"

namespace ipc {

// -- Single collision methods -------------------------------------------------

double SmoothContactPotential::operator()(
    const Collision& collision, const VectorMax12d& positions) const
{
    return collision.weight * collision(positions, params);
}

VectorMax12d SmoothContactPotential::gradient(
    const Collision& collision, const VectorMax12d& positions) const
{
    return collision.weight * collision.gradient(positions, params);
}

MatrixMax12d SmoothContactPotential::hessian(
    const Collision& collision,
    const VectorMax12d& positions,
    const bool project_hessian_to_psd) const
{
    return collision.weight * collision.hessian(positions, params, project_hessian_to_psd);
}

} // namespace ipc