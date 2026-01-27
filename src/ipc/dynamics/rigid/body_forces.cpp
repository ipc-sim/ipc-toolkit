#include "body_forces.hpp"

#include <ipc/dynamics/rigid/rigid_body.hpp>
#include <ipc/geometry/normal.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc::rigid {

void BodyForces::update(const RigidBodies& bodies)
{
    const double dt_sq = time_integrator->dt * time_integrator->dt;

    forces.resize(bodies.num_bodies());
    torques.resize(bodies.num_bodies());

    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, bodies.num_bodies()),
        [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i < r.end(); ++i) {
                forces[i] = -dt_sq
                    * (bodies[i].mass() * gravity()
                       + bodies[i].external_force().position);

                const auto& torque = bodies[i].external_force().rotation;

                // Add external torques to the predicted pose
                if (torque.size() == 3) {
                    const auto& Q = time_integrator->pose(i).rotation;
                    // Transform the world space torque into body space
                    const Eigen::Matrix3d Tau =
                        Q.transpose() * cross_product_matrix(torque);
                    torques[i] = -dt_sq * Tau;
                } else {
                    assert(torque.size() == 1);
                    torques[i].resize(1, 1);
                    torques[i](0, 0) = -dt_sq * torque(0);
                }
            }
        });
}

// ---- Cumulative functions ---------------------------------------------------

double BodyForces::operator()(
    const RigidBodies& bodies, Eigen::ConstRef<Eigen::VectorXd> x)
{
    const int ndof = x.size() / bodies.num_bodies();

    double energy = 0.0;
    for (size_t i = 0; i < bodies.num_bodies(); ++i) {
        energy += operator()(
            bodies[i], x.segment(i * ndof, ndof), forces[i], torques[i]);
    }
    return energy;
}

Eigen::VectorXd BodyForces::gradient(
    const RigidBodies& bodies, Eigen::ConstRef<Eigen::VectorXd> x)
{
    const int ndof = x.size() / bodies.num_bodies();

    Eigen::VectorXd grad = Eigen::VectorXd::Zero(x.size());
    for (size_t i = 0; i < bodies.num_bodies(); ++i) {
        grad.segment(i * ndof, ndof) = gradient(
            bodies[i], x.segment(i * ndof, ndof), forces[i], torques[i]);
    }
    return grad;
}

Eigen::MatrixXd BodyForces::hessian(
    const RigidBodies& bodies,
    Eigen::ConstRef<Eigen::VectorXd> x,
    const PSDProjectionMethod project_hessian_to_psd)
{
    const int ndof = x.size() / bodies.num_bodies();

    Eigen::MatrixXd hess(x.size(), x.size());
    for (size_t i = 0; i < bodies.num_bodies(); ++i) {
        hess.block(i * ndof, i * ndof, ndof, ndof) = hessian(
            bodies[i], x.segment(i * ndof, ndof), forces[i], torques[i],
            project_hessian_to_psd);
    }
    return hess;
}

// ---- Per-body functions -----------------------------------------------------

double BodyForces::operator()(
    const RigidBody& body,
    Eigen::ConstRef<VectorMax6d> x,
    Eigen::ConstRef<VectorMax3d> q_hat,
    Eigen::ConstRef<MatrixMax3d> Q_hat) const
{
    double energy = 0.0;

    // if (!body.is_dof_fixed.head(pose.pos_ndof()).all())
    {
        energy += x.head(forces.size()).dot(q_hat);
    }

    // Rotational energy
    // if (!body.is_dof_fixed.tail(pose.rot_ndof()).all())
    {
        if (Q_hat.size() == 9) {
            const Eigen::Matrix3d Q = rotation_vector_to_matrix(x.tail<3>());
            energy += (Q.transpose() * Q_hat).trace();
        } else {
            assert(Q_hat.size() == 1);
            energy += x(2) * Q_hat(0, 0);
        }
    }

    return energy;
}

VectorMax6d BodyForces::gradient(
    const RigidBody& body,
    Eigen::ConstRef<VectorMax6d> x,
    Eigen::ConstRef<VectorMax3d> q_hat,
    Eigen::ConstRef<MatrixMax3d> Q_hat) const
{
    VectorMax6d grad = VectorMax6d::Zero(x.size());

    // if (!body.is_dof_fixed.head(pose.pos_ndof()).all())
    {
        grad.head(gravity().size()) = q_hat;
    }

    // Rotational energy
    // if (!body.is_dof_fixed.tail(pose.rot_ndof()).all())
    {
        if (Q_hat.size() == 9) {
            const Eigen::Matrix<double, 9, 3> dQ_dx =
                rotation_vector_to_matrix_jacobian(x.tail<3>());
            grad.tail<3>() = dQ_dx.transpose() * Q_hat.reshaped();
        } else {
            assert(Q_hat.size() == 1);
            grad(2) = Q_hat(0, 0);
        }
    }

    return grad;
}

MatrixMax6d BodyForces::hessian(
    const RigidBody& body,
    Eigen::ConstRef<VectorMax6d> x,
    Eigen::ConstRef<VectorMax3d> q_hat,
    Eigen::ConstRef<MatrixMax3d> Q_hat,
    const PSDProjectionMethod project_hessian_to_psd) const
{
    MatrixMax6d hess = MatrixMax6d::Zero(x.size(), x.size());

    // Rotational energy
    // if (!body.is_dof_fixed.tail(pose.rot_ndof()).all())
    {
        if (Q_hat.size() == 9) {
            const Eigen::Matrix<double, 9, 9> d2Q_dx2 =
                rotation_vector_to_matrix_hessian(x.tail<3>());

            const Eigen::Matrix3d hess_rotation =
                (d2Q_dx2.transpose() * Q_hat.reshaped()).reshaped(3, 3);

            hess.bottomRightCorner<3, 3>() =
                project_to_psd(hess_rotation, project_hessian_to_psd);
        }
    }

    return hess;
}

} // namespace ipc::rigid