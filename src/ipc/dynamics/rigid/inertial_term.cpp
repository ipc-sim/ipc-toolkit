#include "inertial_term.hpp"

#include <ipc/dynamics/rigid/rigid_body.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc::rigid {

namespace {
    inline Eigen::Matrix3d
    cross_product_matrix(Eigen::ConstRef<Eigen::Vector3d> x)
    {
        Eigen::Matrix3d X;
        X << 0, -x.z(), x.y(), //
            x.z(), 0, -x.x(),  //
            -x.y(), x.x(), 0;
        return X;
    }
} // namespace

void InertialTerm::update(const RigidBodies& bodies)
{
    // Update the predicted poses based on the current time integrator state
    predicted_poses = time_integrator->predicted_pose();

    // Gravity in the y-direction
    const double dt_sq = time_integrator->dt * time_integrator->dt;

    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, predicted_poses.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i < r.end(); ++i) {
                // Add gravity to the predicted pose
                // TODO: Make this configurable
                predicted_poses[i].position.y() += dt_sq * -9.81;

                const auto& force = bodies[i].external_force().position;
                const auto& torque = bodies[i].external_force().rotation;

                // Add external forces to the predicted pose
                if (!force.isZero()) {
                    predicted_poses[i].position +=
                        dt_sq * force / bodies[i].mass();
                }

                // Add external torques to the predicted pose
                if (!torque.isZero()) {
                    if (torque.size() == 3) {
                        const auto& Q = time_integrator->pose(i).rotation;
                        // Transform the world space torque into body space
                        const Eigen::Matrix3d Tau =
                            Q.transpose() * cross_product_matrix(torque);
                        predicted_poses[i].rotation +=
                            dt_sq * bodies[i].J().inverse() * Tau;
                    } else {
                        assert(torque.size() == 1);
                        predicted_poses[i].rotation(0) += dt_sq * torque(0)
                            / bodies[i].moment_of_inertia()(0);
                    }
                }
            }
        });
}

// ---- Cumulative functions ---------------------------------------------------

double InertialTerm::operator()(
    const RigidBodies& bodies, Eigen::ConstRef<Eigen::VectorXd> x)
{
    const int ndof = x.size() / bodies.num_bodies();

    double energy = 0.0;
    for (size_t i = 0; i < bodies.num_bodies(); ++i) {
        energy += operator()(
            bodies[i], x.segment(i * ndof, ndof), predicted_poses[i].position,
            predicted_poses[i].rotation);
    }
    return energy;
}

Eigen::VectorXd InertialTerm::gradient(
    const RigidBodies& bodies, Eigen::ConstRef<Eigen::VectorXd> x)
{
    const int ndof = x.size() / bodies.num_bodies();

    Eigen::VectorXd grad = Eigen::VectorXd::Zero(x.size());
    for (size_t i = 0; i < bodies.num_bodies(); ++i) {
        grad.segment(i * ndof, ndof) = gradient(
            bodies[i], x.segment(i * ndof, ndof), predicted_poses[i].position,
            predicted_poses[i].rotation);
    }
    return grad;
}

Eigen::MatrixXd InertialTerm::hessian(
    const RigidBodies& bodies, Eigen::ConstRef<Eigen::VectorXd> x)
{
    const int ndof = x.size() / bodies.num_bodies();

    Eigen::MatrixXd hess(x.size(), x.size());
    for (size_t i = 0; i < bodies.num_bodies(); ++i) {
        hess.block(i * ndof, i * ndof, ndof, ndof) = hessian(
            bodies[i], x.segment(i * ndof, ndof), predicted_poses[i].position,
            predicted_poses[i].rotation);
    }
    return hess;
}

// ---- Per-body functions -----------------------------------------------------

double InertialTerm::operator()(
    const RigidBody& body,
    Eigen::ConstRef<VectorMax6d> x,
    Eigen::ConstRef<VectorMax3d> q_hat,
    Eigen::ConstRef<MatrixMax3d> Q_hat) const
{
    double energy = 0.0;

    // Linear inertia
    // if (!body.is_dof_fixed.head(pose.pos_ndof()).all())
    {
        // ½m‖q - q̂‖²
        energy +=
            0.5 * body.mass() * (x.head(q_hat.size()) - q_hat).squaredNorm();
    }

    // Rotational energy
    // if (!body.is_dof_fixed.tail(pose.rot_ndof()).all())
    {
        if (q_hat.size() == 3) {
            const Eigen::Matrix3d Q = rotation_vector_to_matrix(x.tail<3>());

            // ½tr((Q - Q̂) J (Q - Q̂)ᵀ)
            const Eigen::Matrix3d dQ = Q - Q_hat;
            energy += 0.5 * (dQ * body.J() * dQ.transpose()).trace();
        } else {
            assert(q_hat.size() == 2);
            assert(x.size() == 3);
            assert(Q_hat.size() == 1);
            // ½ I (θ-θ̂)²
            const double dtheta = x(2) - Q_hat(0, 0);
            energy += 0.5 * body.moment_of_inertia()(0) * dtheta * dtheta;
        }
    }

    return energy;
}

VectorMax6d InertialTerm::gradient(
    const RigidBody& body,
    Eigen::ConstRef<VectorMax6d> x,
    Eigen::ConstRef<VectorMax3d> q_hat,
    Eigen::ConstRef<MatrixMax3d> Q_hat) const
{
    VectorMax6d grad = VectorMax6d::Zero(x.size());

    // Linear inertia
    // if (!body.is_dof_fixed.head(pose.pos_ndof()).all())
    {
        // m (q - q̂)
        grad.head(q_hat.size()) = body.mass() * (x.head(q_hat.size()) - q_hat);
    }

    // Rotational energy
    // if (!body.is_dof_fixed.tail(pose.rot_ndof()).all())
    {
        if (q_hat.size() == 3) {
            const Eigen::Matrix3d Q = rotation_vector_to_matrix(x.tail<3>());

            // ∂f​/∂xₖ = ∑ᵢ∑ⱼ ∂E/∂Qᵢⱼ ∂Qᵢⱼ/∂xₖ
            const Eigen::Matrix<double, 9, 3> dQ_dx =
                rotation_vector_to_matrix_jacobian(x.tail<3>());
            const Eigen::Vector<double, 9> dE_dQ =
                ((Q - Q_hat) * body.J()).reshaped();
            grad.tail<3>() = dQ_dx.transpose() * dE_dQ;
        } else {
            assert(q_hat.size() == 2);
            assert(x.size() == 3);
            assert(Q_hat.size() == 1);
            // I (θ - θ̂)
            grad(2) = body.moment_of_inertia()(0) * (x(2) - Q_hat(0, 0));
        }
    }

    return grad;
}

MatrixMax6d InertialTerm::hessian(
    const RigidBody& body,
    Eigen::ConstRef<VectorMax6d> x,
    Eigen::ConstRef<VectorMax3d> q_hat,
    Eigen::ConstRef<MatrixMax3d> Q_hat) const
{
    MatrixMax6d hess = MatrixMax6d::Zero(x.size(), x.size());

    // Linear inertia
    // if (!body.is_dof_fixed.head(pose.pos_ndof()).all())
    {
        // m (q - q̂)
        hess.topLeftCorner(q_hat.size(), q_hat.size()).diagonal().array() =
            body.mass();
    }

    // Rotational energy
    // if (!body.is_dof_fixed.tail(pose.rot_ndof()).all())
    {
        if (q_hat.size() == 3) {
            const Eigen::Matrix3d Q = rotation_vector_to_matrix(x.tail<3>());
            const Eigen::Matrix<double, 9, 3> dQ_dx =
                rotation_vector_to_matrix_jacobian(x.tail<3>());
            const Eigen::Matrix<double, 9, 9> d2Q_dx2 =
                rotation_vector_to_matrix_hessian(x.tail<3>());

            const Eigen::Vector<double, 9> dE_dQ =
                ((Q - Q_hat) * body.J()).reshaped();

            // ∂²E/∂Q² = J ⊗ I
            Eigen::Matrix<double, 9, 9> d2E_dQ2 =
                Eigen::Matrix<double, 9, 9>::Zero();
            d2E_dQ2.diagonal().segment<3>(0).array() = body.J()(0, 0);
            d2E_dQ2.diagonal().segment<3>(3).array() = body.J()(1, 1);
            d2E_dQ2.diagonal().segment<3>(6).array() = body.J()(2, 2);

            // (3x3) = (3×9)(9×9)(9x3) + mat((9x9)(9x1))
            hess.bottomRightCorner<3, 3>() = dQ_dx.transpose() * d2E_dQ2 * dQ_dx
                + (d2Q_dx2.transpose() * dE_dQ).reshaped(3, 3);
        } else {
            assert(q_hat.size() == 2);
            assert(x.size() == 3);
            assert(Q_hat.size() == 1);
            // I (θ - θ̂)
            hess(2, 2) = body.moment_of_inertia()(0);
        }
    }

    return hess;
}

} // namespace ipc::rigid