#include "augmented_lagrangian.hpp"

#include <ipc/utils/logger.hpp>

namespace ipc::affine {

namespace {
    /// W = J ⊗ I: diagonal of the angular ABD mass block for A(k, j) DOFs
    /// (DOF layout ndof·i + dim + k + dim·j ⇒ weight J(j)).
    Eigen::VectorXd angular_weights(const rigid::RigidBody& body) // NOLINT
    {
        const int dim = body.R0().rows();
        Eigen::VectorXd w(dim * dim);
        for (int j = 0; j < dim; ++j) {
            w.segment(dim * j, dim).array() = body.J().diagonal()(j);
        }
        return w;
    }

    /// Progress η = 1 − √(d²/d₀²) with an absolute floor (see the rigid AL).
    double progress(
        const double distance_sq,
        const double start_distance_sq,
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
    const rigid::RigidBodies& bodies,
    Eigen::ConstRef<Eigen::VectorXd> x0,
    const std::vector<rigid::Pose>& targets)
{
    assert(targets.size() == bodies.num_bodies());
    m_dim = bodies.dim();
    assert(
        x0.size()
        == (m_dim + m_dim * m_dim) * Eigen::Index(bodies.num_bodies()));

    m_x0 = x0;
    m_targets = targets;

    m_target_rotations.resize(targets.size());
    m_lambda.resize(targets.size());
    m_lambda_A.resize(targets.size());
    m_num_kinematic = 0;
    for (size_t i = 0; i < targets.size(); ++i) {
        m_lambda[i] = VectorMax3d::Zero(m_dim);
        m_lambda_A[i] = Eigen::VectorXd::Zero(m_dim * m_dim);
        // Only kinematic bodies have a meaningful target rotation; the rest
        // are left default (never read — see the KINEMATIC guards below).
        if (bodies[i].type() == rigid::RigidBody::Type::KINEMATIC) {
            m_target_rotations[i] = targets[i].rotation_matrix();
            ++m_num_kinematic;
        }
    }

    m_kappa_q = m_params.initial_penalty;
    m_kappa_A = m_params.initial_penalty;

    m_linear_satisfied = m_num_kinematic == 0;
    m_angular_satisfied = m_num_kinematic == 0;
}

VectorMax12d AugmentedLagrangian::target_dof(const size_t body_id) const
{
    VectorMax12d y(m_dim + m_dim * m_dim);
    y.head(m_dim) = m_targets[body_id].position;
    y.tail(m_dim * m_dim) = m_target_rotations[body_id].reshaped();
    return y;
}

// ---- Progress metrics ---------------------------------------------------

double AugmentedLagrangian::linear_progress(
    const rigid::RigidBodies& bodies, Eigen::ConstRef<Eigen::VectorXd> x) const
{
    double distance_sq = 0, start_distance_sq = 0;
    for (size_t i = 0; i < bodies.num_bodies(); ++i) {
        if (bodies[i].type() != rigid::RigidBody::Type::KINEMATIC) {
            continue;
        }
        const int ndof = m_dim + m_dim * m_dim;
        const auto& q_hat = m_targets[i].position;
        distance_sq += (q_hat - x.segment(ndof * i, m_dim)).squaredNorm();
        start_distance_sq +=
            (q_hat - m_x0.segment(ndof * i, m_dim)).squaredNorm();
    }

    double scale = 0;
    for (size_t i = 0; i < bodies.num_bodies(); ++i) {
        if (bodies[i].type() == rigid::RigidBody::Type::KINEMATIC) {
            scale = std::max(scale, bodies[i].bounding_radius());
        }
    }
    return progress(
        distance_sq, start_distance_sq,
        (1.0 - m_params.satisfied_progress) * std::max(scale, 1e-12));
}

double AugmentedLagrangian::angular_progress(
    const rigid::RigidBodies& bodies, Eigen::ConstRef<Eigen::VectorXd> x) const
{
    double distance_sq = 0, start_distance_sq = 0;
    for (size_t i = 0; i < bodies.num_bodies(); ++i) {
        if (bodies[i].type() != rigid::RigidBody::Type::KINEMATIC) {
            continue;
        }
        const int ndof = m_dim + m_dim * m_dim;
        const Eigen::VectorXd a_hat = m_target_rotations[i].reshaped();
        distance_sq +=
            (a_hat - x.segment(ndof * i + m_dim, m_dim * m_dim)).squaredNorm();
        start_distance_sq +=
            (a_hat - m_x0.segment(ndof * i + m_dim, m_dim * m_dim))
                .squaredNorm();
    }

    // A matrices have O(1) Frobenius scale.
    return progress(
        distance_sq, start_distance_sq, 1.0 - m_params.satisfied_progress);
}

// ---- Update policy --------------------------------------------------------

void AugmentedLagrangian::update(
    const rigid::RigidBodies& bodies, Eigen::ConstRef<Eigen::VectorXd> x)
{
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
                if (bodies[i].type() != rigid::RigidBody::Type::KINEMATIC) {
                    continue;
                }
                const int ndof = m_dim + m_dim * m_dim;
                m_lambda[i] -= m_kappa_q * std::sqrt(bodies[i].mass())
                    * (x.segment(ndof * i, m_dim) - m_targets[i].position);
            }
        }
        logger().debug(
            "augmented lagrangian (linear): η_q={:g} κ_q={:g} satisfied={}",
            eta_q, m_kappa_q, m_linear_satisfied);
    }

    if (!m_angular_satisfied) {
        const double eta_A = angular_progress(bodies, x); // NOLINT
        if (eta_A >= m_params.satisfied_progress) {
            m_angular_satisfied = true;
        } else if (
            eta_A < m_params.stall_progress
            && m_kappa_A < m_params.max_penalty) {
            m_kappa_A *= 2;
        } else {
            for (size_t i = 0; i < bodies.num_bodies(); ++i) {
                if (bodies[i].type() != rigid::RigidBody::Type::KINEMATIC) {
                    continue;
                }
                // λ_A ← λ_A − κ_A W^{1/2} (a − â)
                const int ndof = m_dim + m_dim * m_dim;
                m_lambda_A[i] -= m_kappa_A
                    * angular_weights(bodies[i])
                          .array()
                          .sqrt()
                          .matrix()
                          .asDiagonal()
                    * (x.segment(ndof * i + m_dim, m_dim * m_dim)
                       - Eigen::VectorXd(m_target_rotations[i].reshaped()));
            }
        }
        logger().debug(
            "augmented lagrangian (angular): η_A={:g} κ_A={:g} satisfied={}",
            eta_A, m_kappa_A, m_angular_satisfied);
    }
}

// ---- Cumulative energy --------------------------------------------------

double AugmentedLagrangian::operator()(
    const rigid::RigidBodies& bodies, Eigen::ConstRef<Eigen::VectorXd> x) const
{
    if (!active()) {
        return 0.0;
    }

    double energy = 0.0;
    for (size_t i = 0; i < bodies.num_bodies(); ++i) {
        if (bodies[i].type() != rigid::RigidBody::Type::KINEMATIC) {
            continue;
        }
        const int ndof = m_dim + m_dim * m_dim;
        const VectorMax3d dq =
            x.segment(ndof * i, m_dim) - m_targets[i].position;
        energy += 0.5 * m_kappa_q * bodies[i].mass() * dq.squaredNorm()
            - std::sqrt(bodies[i].mass()) * m_lambda[i].dot(dq);

        const Eigen::VectorXd da = x.segment(ndof * i + m_dim, m_dim * m_dim)
            - Eigen::VectorXd(m_target_rotations[i].reshaped());
        const Eigen::VectorXd w = angular_weights(bodies[i]);
        energy += 0.5 * m_kappa_A * da.dot(w.cwiseProduct(da))
            - m_lambda_A[i].dot(w.array().sqrt().matrix().cwiseProduct(da));
    }
    return energy;
}

Eigen::VectorXd AugmentedLagrangian::gradient(
    const rigid::RigidBodies& bodies, Eigen::ConstRef<Eigen::VectorXd> x) const
{
    Eigen::VectorXd grad = Eigen::VectorXd::Zero(x.size());
    if (!active()) {
        return grad;
    }

    for (size_t i = 0; i < bodies.num_bodies(); ++i) {
        if (bodies[i].type() != rigid::RigidBody::Type::KINEMATIC) {
            continue;
        }
        const int ndof = m_dim + m_dim * m_dim;
        const VectorMax3d dq =
            x.segment(ndof * i, m_dim) - m_targets[i].position;
        grad.segment(ndof * i, m_dim) = m_kappa_q * bodies[i].mass() * dq
            - std::sqrt(bodies[i].mass()) * m_lambda[i];

        const Eigen::VectorXd da = x.segment(ndof * i + m_dim, m_dim * m_dim)
            - Eigen::VectorXd(m_target_rotations[i].reshaped());
        const Eigen::VectorXd w = angular_weights(bodies[i]);
        grad.segment(ndof * i + m_dim, m_dim * m_dim) =
            m_kappa_A * w.cwiseProduct(da)
            - w.array().sqrt().matrix().cwiseProduct(m_lambda_A[i]);
    }
    return grad;
}

Eigen::SparseMatrix<double> AugmentedLagrangian::hessian(
    const rigid::RigidBodies& bodies, Eigen::ConstRef<Eigen::VectorXd> x) const
{
    Eigen::SparseMatrix<double> hess(x.size(), x.size());
    if (!active()) {
        return hess;
    }

    std::vector<Eigen::Triplet<double>> triplets;
    for (size_t i = 0; i < bodies.num_bodies(); ++i) {
        if (bodies[i].type() != rigid::RigidBody::Type::KINEMATIC) {
            continue;
        }
        const int ndof = m_dim + m_dim * m_dim;
        for (int k = 0; k < m_dim; ++k) {
            triplets.emplace_back(
                ndof * i + k, ndof * i + k, m_kappa_q * bodies[i].mass());
        }
        const Eigen::VectorXd w = angular_weights(bodies[i]);
        for (int k = 0; k < m_dim * m_dim; ++k) {
            triplets.emplace_back(
                ndof * i + m_dim + k, ndof * i + m_dim + k, m_kappa_A * w(k));
        }
    }
    hess.setFromTriplets(triplets.begin(), triplets.end());
    return hess;
}

} // namespace ipc::affine
