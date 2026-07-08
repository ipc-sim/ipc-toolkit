#include "augmented_lagrangian.hpp"

#include <ipc/utils/logger.hpp>

namespace ipc::rigid {

namespace {
    /// J^{1/2}: elementwise sqrt of the diagonal second-moment matrix.
    VectorMax3d sqrt_J(const RigidBody& body) // NOLINT
    {
        return body.J().diagonal().array().sqrt();
    }

    /// Progress η = 1 − √(d²/d₀²) with an absolute floor: when the intended
    /// motion d₀ is below tol (e.g., a translating body with a fixed target
    /// rotation), the relative metric is meaningless — any drift (contact
    /// forces perturb unpinned DOFs) would read as η = -∞. Satisfied iff the
    /// deviation itself is below tol.
    double progress(
        const double distance_sq, const double start_distance_sq,
        const double tol)
    {
        const double tol_sq = tol * tol;
        if (start_distance_sq <= tol_sq) {
            return distance_sq <= tol_sq ? 1.0 : 0.0;
        }
        return 1.0 - std::sqrt(distance_sq / start_distance_sq);
    }
} // namespace

void AugmentedLagrangian::init(
    const RigidBodies& bodies,
    Eigen::ConstRef<Eigen::VectorXd> x0,
    const std::vector<Pose>& targets)
{
    assert(targets.size() == bodies.num_bodies());

    m_start = Pose::to_poses(x0, bodies.dim());
    m_targets = targets;

    m_target_rotations.resize(targets.size());
    m_lambda.resize(targets.size());
    m_Lambda.resize(targets.size());
    m_num_kinematic = 0;
    for (size_t i = 0; i < targets.size(); ++i) {
        m_target_rotations[i] = targets[i].rotation_matrix();
        m_lambda[i] = VectorMax3d::Zero(targets[i].position.size());
        m_Lambda[i] = MatrixMax3d::Zero(
            m_target_rotations[i].rows(), m_target_rotations[i].cols());
        if (bodies[i].type() == RigidBody::Type::KINEMATIC) {
            ++m_num_kinematic;
        }
    }

    m_kappa_q = m_params.initial_penalty;
    m_kappa_Q = m_params.initial_penalty;

    // No kinematic bodies ⇒ nothing to satisfy.
    m_linear_satisfied = m_num_kinematic == 0;
    m_angular_satisfied = m_num_kinematic == 0;
}

// ---- Progress metrics ---------------------------------------------------

double AugmentedLagrangian::linear_progress(
    const RigidBodies& bodies, Eigen::ConstRef<Eigen::VectorXd> x) const
{
    const std::vector<Pose> poses = Pose::to_poses(x, bodies.dim());

    double distance_sq = 0, start_distance_sq = 0;
    for (size_t i = 0; i < bodies.num_bodies(); ++i) {
        if (bodies[i].type() != RigidBody::Type::KINEMATIC) {
            continue;
        }
        distance_sq +=
            (m_targets[i].position - poses[i].position).squaredNorm();
        start_distance_sq +=
            (m_targets[i].position - m_start[i].position).squaredNorm();
    }

    // Length scale for the absolute tolerance: the kinematic bodies' size.
    double scale = 0;
    for (size_t i = 0; i < bodies.num_bodies(); ++i) {
        if (bodies[i].type() == RigidBody::Type::KINEMATIC) {
            scale = std::max(scale, bodies[i].bounding_radius());
        }
    }
    return progress(
        distance_sq, start_distance_sq,
        (1.0 - m_params.satisfied_progress) * std::max(scale, 1e-12));
}

double AugmentedLagrangian::angular_progress(
    const RigidBodies& bodies, Eigen::ConstRef<Eigen::VectorXd> x) const
{
    const std::vector<Pose> poses = Pose::to_poses(x, bodies.dim());

    double distance_sq = 0, start_distance_sq = 0;
    for (size_t i = 0; i < bodies.num_bodies(); ++i) {
        if (bodies[i].type() != RigidBody::Type::KINEMATIC) {
            continue;
        }
        distance_sq +=
            (m_target_rotations[i] - poses[i].rotation_matrix()).squaredNorm();
        start_distance_sq +=
            (m_target_rotations[i] - m_start[i].rotation_matrix())
                .squaredNorm();
    }

    // Rotation matrices have O(1) Frobenius scale.
    return progress(
        distance_sq, start_distance_sq, 1.0 - m_params.satisfied_progress);
}

// ---- Update policy --------------------------------------------------------

void AugmentedLagrangian::update(
    const RigidBodies& bodies, Eigen::ConstRef<Eigen::VectorXd> x)
{
    const std::vector<Pose> poses = Pose::to_poses(x, bodies.dim());

    if (!m_linear_satisfied) {
        const double eta_q = linear_progress(bodies, x);
        if (eta_q >= m_params.satisfied_progress) {
            m_linear_satisfied = true;
        } else if (
            eta_q < m_params.stall_progress
            && m_kappa_q < m_params.max_penalty) {
            m_kappa_q *= 2;
        } else {
            for (size_t i = 0; i < bodies.num_bodies(); ++i) {
                if (bodies[i].type() != RigidBody::Type::KINEMATIC) {
                    continue;
                }
                m_lambda[i] -= m_kappa_q * std::sqrt(bodies[i].mass())
                    * (poses[i].position - m_targets[i].position);
            }
        }
        logger().debug(
            "augmented lagrangian (linear): η_q={:g} κ_q={:g} satisfied={}",
            eta_q, m_kappa_q, m_linear_satisfied);
    }

    if (!m_angular_satisfied) {
        const double eta_Q = angular_progress(bodies, x); // NOLINT
        if (eta_Q >= m_params.satisfied_progress) {
            m_angular_satisfied = true;
        } else if (
            eta_Q < m_params.stall_progress
            && m_kappa_Q < m_params.max_penalty) {
            m_kappa_Q *= 2;
        } else {
            for (size_t i = 0; i < bodies.num_bodies(); ++i) {
                if (bodies[i].type() != RigidBody::Type::KINEMATIC) {
                    continue;
                }
                // Λ ← Λ − κ_Q (Q − Q̂) J^{1/2}
                m_Lambda[i] -= m_kappa_Q
                    * (poses[i].rotation_matrix() - m_target_rotations[i])
                    * sqrt_J(bodies[i]).asDiagonal();
            }
        }
        logger().debug(
            "augmented lagrangian (angular): η_Q={:g} κ_Q={:g} satisfied={}",
            eta_Q, m_kappa_Q, m_angular_satisfied);
    }
}

// ---- Per-body energy --------------------------------------------------------

double AugmentedLagrangian::operator()(
    const size_t body_id,
    const RigidBody& body,
    Eigen::ConstRef<VectorMax6d> x) const
{
    if (!active() || body.type() != RigidBody::Type::KINEMATIC) {
        return 0.0;
    }

    const int dim = m_targets[body_id].position.size();

    double energy = 0.0;

    // Linear: (κ_q/2) m ‖q − q̂‖² − √m λᵀ(q − q̂)
    const VectorMax3d dq = x.head(dim) - m_targets[body_id].position;
    energy += 0.5 * m_kappa_q * body.mass() * dq.squaredNorm()
        - std::sqrt(body.mass()) * m_lambda[body_id].dot(dq);

    // Angular: (κ_Q/2) tr((Q − Q̂) J (Q − Q̂)ᵀ) − tr(Λᵀ (Q − Q̂) J^{1/2})
    if (dim == 3) {
        const Eigen::Matrix3d dQ = rotation_vector_to_matrix(x.tail<3>())
            - m_target_rotations[body_id];
        energy += 0.5 * m_kappa_Q
                * (dQ * body.J() * dQ.transpose()).trace()
            - (m_Lambda[body_id].transpose() * dQ
               * sqrt_J(body).asDiagonal())
                  .trace();
    } else {
        // 2D: (κ_Q/2) I (θ − θ̂)² − √I Λ (θ − θ̂)
        const double dtheta = x(2) - m_targets[body_id].rotation(0);
        const double I = body.moment_of_inertia()(0); // NOLINT
        energy += 0.5 * m_kappa_Q * I * dtheta * dtheta
            - std::sqrt(I) * m_Lambda[body_id](0, 0) * dtheta;
    }

    return energy;
}

VectorMax6d AugmentedLagrangian::gradient(
    const size_t body_id,
    const RigidBody& body,
    Eigen::ConstRef<VectorMax6d> x) const
{
    if (!active() || body.type() != RigidBody::Type::KINEMATIC) {
        return VectorMax6d::Zero(x.size());
    }

    const int dim = m_targets[body_id].position.size();

    VectorMax6d grad = VectorMax6d::Zero(x.size());

    // ∇_q = κ_q m (q − q̂) − √m λ
    const VectorMax3d dq = x.head(dim) - m_targets[body_id].position;
    grad.head(dim) = m_kappa_q * body.mass() * dq
        - std::sqrt(body.mass()) * m_lambda[body_id];

    if (dim == 3) {
        // Matrix-space gradient G = κ_Q (Q − Q̂) J − Λ J^{1/2}, then chain
        // rule through dQ/dθ (as in InertialTerm).
        const Eigen::Matrix3d dQ = rotation_vector_to_matrix(x.tail<3>())
            - m_target_rotations[body_id];
        const Eigen::Matrix3d G = m_kappa_Q * dQ * body.J()
            - m_Lambda[body_id] * sqrt_J(body).asDiagonal();

        const Eigen::Matrix<double, 9, 3> dQ_dx =
            rotation_vector_to_matrix_jacobian(x.tail<3>());
        grad.tail<3>() = dQ_dx.transpose() * G.reshaped();
    } else {
        const double dtheta = x(2) - m_targets[body_id].rotation(0);
        const double I = body.moment_of_inertia()(0); // NOLINT
        grad(2) =
            m_kappa_Q * I * dtheta - std::sqrt(I) * m_Lambda[body_id](0, 0);
    }

    return grad;
}

MatrixMax6d AugmentedLagrangian::hessian(
    const size_t body_id,
    const RigidBody& body,
    Eigen::ConstRef<VectorMax6d> x,
    const PSDProjectionMethod project_hessian_to_psd) const
{
    if (!active() || body.type() != RigidBody::Type::KINEMATIC) {
        return MatrixMax6d::Zero(x.size(), x.size());
    }

    const int dim = m_targets[body_id].position.size();

    MatrixMax6d hess = MatrixMax6d::Zero(x.size(), x.size());

    // ∇²_q = κ_q m I (constant, PSD)
    hess.topLeftCorner(dim, dim).diagonal().array() =
        m_kappa_q * body.mass();

    if (dim == 3) {
        const Eigen::Matrix3d dQ = rotation_vector_to_matrix(x.tail<3>())
            - m_target_rotations[body_id];
        const Eigen::Matrix3d G = m_kappa_Q * dQ * body.J()
            - m_Lambda[body_id] * sqrt_J(body).asDiagonal();

        const Eigen::Matrix<double, 9, 3> dQ_dx =
            rotation_vector_to_matrix_jacobian(x.tail<3>());
        const Eigen::Matrix<double, 9, 9> d2Q_dx2 =
            rotation_vector_to_matrix_hessian(x.tail<3>());

        // ∂²E/∂Q² = κ_Q (J ⊗ I₃) (the Λ term is linear in Q)
        Eigen::Matrix<double, 9, 9> d2E_dQ2 =
            Eigen::Matrix<double, 9, 9>::Zero();
        d2E_dQ2.diagonal().segment<3>(0).array() =
            m_kappa_Q * body.J().diagonal()(0);
        d2E_dQ2.diagonal().segment<3>(3).array() =
            m_kappa_Q * body.J().diagonal()(1);
        d2E_dQ2.diagonal().segment<3>(6).array() =
            m_kappa_Q * body.J().diagonal()(2);

        const Eigen::Matrix3d hess_rotation =
            dQ_dx.transpose() * d2E_dQ2 * dQ_dx
            + (d2Q_dx2.transpose() * G.reshaped()).reshaped(3, 3);

        hess.bottomRightCorner<3, 3>() =
            project_to_psd(hess_rotation, project_hessian_to_psd);
    } else {
        hess(2, 2) = m_kappa_Q * body.moment_of_inertia()(0);
    }

    return hess;
}

} // namespace ipc::rigid
