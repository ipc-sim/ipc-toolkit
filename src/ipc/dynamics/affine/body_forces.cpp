#include "body_forces.hpp"

#include <ipc/dynamics/time_integration/time_integrator.hpp>
#include <ipc/geometry/normal.hpp>

namespace ipc::affine {

BodyForces::BodyForces(
    const std::shared_ptr<const dynamics::ImplicitTimeIntegrator>&
        _time_integrator)
    : time_integrator(_time_integrator)
{
}

void BodyForces::update(const rigid::RigidBodies& bodies)
{
    assert(bodies.dim() == 3); // TODO: Support 2D affine bodies

    const double dt_sq = time_integrator->acceleration_scaling();

    m_wrench = Eigen::VectorXd::Zero(12 * bodies.num_bodies());
    for (size_t i = 0; i < bodies.num_bodies(); ++i) {
        m_wrench.segment<3>(12 * i) = -dt_sq
            * (bodies[i].mass() * gravity()
               + bodies[i].external_force().position);

        const auto& torque = bodies[i].external_force().rotation;
        if (!torque.isZero()) {
            // Transform the world space torque into body space (as in
            // rigid::BodyForces::update)
            const auto& A = time_integrator->pose(i).rotation;
            const Eigen::Matrix3d tau =
                A.transpose() * cross_product_matrix(torque);
            m_wrench.segment<9>(12 * i + 3) = -dt_sq * tau.reshaped();
        }
    }
}

} // namespace ipc::affine
