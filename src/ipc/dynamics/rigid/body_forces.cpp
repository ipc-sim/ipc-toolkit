#include "body_forces.hpp"

#include <ipc/dynamics/rigid/time_integrator.hpp>
#include <ipc/geometry/normal.hpp>

namespace ipc::rigid {

void BodyForces::update(const RigidBodies& bodies)
{
    const double dt_sq = time_integrator->dt * time_integrator->dt;

    m_forces.resize(bodies.num_bodies());
    m_torques.resize(bodies.num_bodies());

    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, bodies.num_bodies()),
        [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i < r.end(); ++i) {
                m_forces[i] = -dt_sq
                    * (bodies[i].mass() * gravity()
                       + bodies[i].external_force().position);

                const auto& torque = bodies[i].external_force().rotation;

                // Add external torques to the predicted pose
                if (torque.size() == 3) {
                    const auto& Q = time_integrator->pose(i).rotation;
                    // Transform the world space torque into body space
                    const Eigen::Matrix3d Tau =
                        Q.transpose() * cross_product_matrix(torque);
                    m_torques[i] = -dt_sq * Tau;
                } else {
                    assert(torque.size() == 1);
                    m_torques[i].resize(1, 1);
                    m_torques[i](0, 0) = -dt_sq * torque(0);
                }
            }
        });
}

// ---- Per-body functions -----------------------------------------------------

double BodyForces::operator()(
    const size_t body_id,
    const RigidBody& body,
    Eigen::ConstRef<VectorMax6d> x) const
{
    const auto& force = forces()[body_id];
    const auto& torque = torques()[body_id];

    double energy = 0.0;

    // if (!body.is_dof_fixed.head(pose.pos_ndof()).all())
    if (!force.isZero()) {
        energy += x.head(force.size()).dot(force);
    }

    // Rotational energy
    // if (!body.is_dof_fixed.tail(pose.rot_ndof()).all())
    if (!torque.isZero()) {
        if (torque.size() == 9) {
            const Eigen::Matrix3d Q = rotation_vector_to_matrix(x.tail<3>());
            energy += (Q.transpose() * torque).trace();
        } else {
            assert(torque.size() == 1);
            energy += x(2) * torque(0, 0);
        }
    }

    return energy;
}

VectorMax6d BodyForces::gradient(
    const size_t body_id,
    const RigidBody& body,
    Eigen::ConstRef<VectorMax6d> x) const
{
    const auto& force = forces()[body_id];
    const auto& torque = torques()[body_id];

    VectorMax6d grad = VectorMax6d::Zero(x.size());

    // if (!body.is_dof_fixed.head(pose.pos_ndof()).all())
    if (!force.isZero()) {
        grad.head(gravity().size()) = force;
    }

    // Rotational energy
    // if (!body.is_dof_fixed.tail(pose.rot_ndof()).all())
    if (!torque.isZero()) {
        if (torque.size() == 9) {
            const Eigen::Matrix<double, 9, 3> dQ_dx =
                rotation_vector_to_matrix_jacobian(x.tail<3>());
            grad.tail<3>() = dQ_dx.transpose() * torque.reshaped();
        } else {
            assert(torque.size() == 1);
            grad(2) = torque(0, 0);
        }
    }

    return grad;
}

MatrixMax6d BodyForces::hessian(
    const size_t body_id,
    const RigidBody& body,
    Eigen::ConstRef<VectorMax6d> x,
    const PSDProjectionMethod project_hessian_to_psd) const
{
    const auto& torque = torques()[body_id];

    MatrixMax6d hess = MatrixMax6d::Zero(x.size(), x.size());

    // Rotational energy
    // if (!body.is_dof_fixed.tail(pose.rot_ndof()).all())
    {
        if (torque.size() == 9 && !torque.isZero()) {
            const Eigen::Matrix<double, 9, 9> d2Q_dx2 =
                rotation_vector_to_matrix_hessian(x.tail<3>());

            const Eigen::Matrix3d hess_rotation =
                (d2Q_dx2.transpose() * torque.reshaped()).reshaped(3, 3);

            hess.bottomRightCorner<3, 3>() =
                project_to_psd(hess_rotation, project_hessian_to_psd);
        }
    }

    return hess;
}

} // namespace ipc::rigid