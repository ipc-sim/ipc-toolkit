#include "body_forces.hpp"

#include <ipc/geometry/normal.hpp>

namespace ipc::affine {

void BodyForces::update(
    const rigid::RigidBodies& bodies,
    const std::vector<affine::Pose>& poses)
{
    assert(poses.size() == bodies.num_bodies());

    const int dim = bodies.dim();
    const int ndof = dim + dim * dim;

    m_wrench = Eigen::VectorXd::Zero(ndof * bodies.num_bodies());
    for (size_t i = 0; i < bodies.num_bodies(); ++i) {
        if (!bodies[i].is_dynamic()) {
            continue; // no gravity or external forces
        }
        m_wrench.segment(ndof * i, dim) =
            -(bodies[i].mass() * gravity().head(dim)
              + bodies[i].external_force().position);

        const auto& torque = bodies[i].external_force().rotation;
        if (!torque.isZero()) {
            // Transform the world space torque into body space (as in
            // rigid::BodyForces::update)
            const auto& A = poses[i].rotation;
            if (dim == 3) {
                const Eigen::Matrix3d tau =
                    A.transpose() * cross_product_matrix(torque);
                m_wrench.segment<9>(ndof * i + 3) = -tau.reshaped();
            } else {
                // 2D: the skew generator is S = [[0, -1], [1, 0]].
                Eigen::Matrix2d skew; // NOLINT
                skew << 0, -torque(0), torque(0), 0;
                const Eigen::Matrix2d tau = A.transpose() * skew;
                m_wrench.segment<4>(ndof * i + 2) = -tau.reshaped();
            }
        }
    }
}

} // namespace ipc::affine
