#include "body_forces.hpp"

#include <ipc/geometry/normal.hpp>

namespace ipc::affine {

void BodyForces::update(
    const rigid::RigidBodies& bodies, const std::vector<affine::Pose>& poses)
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
        if (torque.isZero()) {
            continue;
        }

        // We pull the world-space torque τ back to the affine matrix A with
        // the linear energy E(A) = −⟨½ [τ]× A₀, A⟩, where A₀ is the step-start
        // A. Perturbing A₀ by a small world rotation δφ, A = (I + [δφ]×) A₀,
        // gives δE = −tr(½ A₀ᵀ [τ]×ᵀ [δφ]× A₀) = −τ·δφ, so the generalized
        // torque is exactly τ regardless of the body's orientation. The ½
        // comes from tr([a]×ᵀ [b]×) = 2 a·b.
        const auto& A = poses[i].rotation;
        if (dim == 3) {
            const Eigen::Matrix3d W =
                0.5 * cross_product_matrix(torque) * Eigen::Matrix3d(A);
            m_wrench.segment<9>(ndof * i + 3) = -W.reshaped();
        } else {
            const Eigen::Matrix2d W =
                0.5 * cross_product_matrix(torque(0)) * Eigen::Matrix2d(A);
            m_wrench.segment<4>(ndof * i + 2) = -W.reshaped();
        }
    }
}

} // namespace ipc::affine
