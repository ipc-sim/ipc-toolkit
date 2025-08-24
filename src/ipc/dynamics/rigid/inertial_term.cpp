#include "inertial_term.hpp"

#include <ipc/dynamics/rigid/rigid_body.hpp>
#include <ipc/utils/eigen_ext.hpp>

#include <iostream>

namespace ipc::rigid {

void InertialTerm::update(const RigidBodies& bodies)
{
    // Update the predicted poses based on the current time integrator state
    m_predicted_poses = time_integrator->predicted_pose();
}

// ---- Cumulative functions ---------------------------------------------------

double InertialTerm::operator()(
    const RigidBodies& bodies, Eigen::ConstRef<Eigen::VectorXd> x)
{
    assert(predicted_poses().size() == bodies.num_bodies());

    const int ndof = x.size() / bodies.num_bodies();

    double energy = 0.0;
    for (size_t i = 0; i < bodies.num_bodies(); ++i) {
        energy += operator()(
            bodies[i], x.segment(i * ndof, ndof), predicted_poses()[i].position,
            predicted_poses()[i].rotation);
    }
    return energy;
}

Eigen::VectorXd InertialTerm::gradient(
    const RigidBodies& bodies, Eigen::ConstRef<Eigen::VectorXd> x)
{
    assert(predicted_poses().size() == bodies.num_bodies());

    const int ndof = x.size() / bodies.num_bodies();

    Eigen::VectorXd grad = Eigen::VectorXd::Zero(x.size());
    for (size_t i = 0; i < bodies.num_bodies(); ++i) {
        grad.segment(i * ndof, ndof) = gradient(
            bodies[i], x.segment(i * ndof, ndof), predicted_poses()[i].position,
            predicted_poses()[i].rotation);
    }
    return grad;
}

Eigen::MatrixXd InertialTerm::hessian(
    const RigidBodies& bodies,
    Eigen::ConstRef<Eigen::VectorXd> x,
    const PSDProjectionMethod project_hessian_to_psd)
{
    assert(predicted_poses().size() == bodies.num_bodies());

    const int ndof = x.size() / bodies.num_bodies();

    Eigen::MatrixXd hess(x.size(), x.size());
    for (size_t i = 0; i < bodies.num_bodies(); ++i) {
        hess.block(i * ndof, i * ndof, ndof, ndof) = hessian(
            bodies[i], x.segment(i * ndof, ndof), predicted_poses()[i].position,
            predicted_poses()[i].rotation, project_hessian_to_psd);
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
    Eigen::ConstRef<MatrixMax3d> Q_hat,
    const PSDProjectionMethod project_hessian_to_psd) const
{
    MatrixMax6d hess = MatrixMax6d::Zero(x.size(), x.size());

    // Linear inertia
    // if (!body.is_dof_fixed.head(pose.pos_ndof()).all())
    {
        // m (q - q̂)
        hess.topLeftCorner(q_hat.size(), q_hat.size()).diagonal().array() =
            body.mass();

        // NOTE: No need to project to PSD here, as the Hessian is already
        // positive semi-definite (diagonal matrix with positive entries).
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
            d2E_dQ2.diagonal().segment<3>(0).array() = body.J().diagonal()(0);
            d2E_dQ2.diagonal().segment<3>(3).array() = body.J().diagonal()(1);
            d2E_dQ2.diagonal().segment<3>(6).array() = body.J().diagonal()(2);

            // (3x3) = (3×9)(9×9)(9x3) + mat((9x9)(9x1))
            Eigen::Matrix3d hess_rotation = dQ_dx.transpose() * d2E_dQ2 * dQ_dx
                + (d2Q_dx2.transpose() * dE_dQ).reshaped(3, 3);

            hess.bottomRightCorner<3, 3>() =
                project_to_psd(hess_rotation, project_hessian_to_psd);

        } else {
            assert(q_hat.size() == 2);
            assert(x.size() == 3);
            assert(Q_hat.size() == 1);
            // I (θ - θ̂)
            hess(2, 2) = body.moment_of_inertia()(0);

            // NOTE: No need to project to PSD here, as the Hessian is already
            // positive semi-definite (single positive value).
            assert(body.moment_of_inertia()(0) > 0.0);
        }
    }

    return hess;
}

} // namespace ipc::rigid