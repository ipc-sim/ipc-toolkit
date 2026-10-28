#include "augmented_lagrangian.hpp"

#include <ipc/utils/logger.hpp>

namespace ipc::affine {

namespace {
    /// W = J ⊗ I: diagonal of the angular ABD mass block for A(k, j) DOFs
    /// (DOF layout ndof·i + dim + k + dim·j ⇒ weight J(j)).
    Eigen::VectorXd angular_weights(const rigid::RigidBody& body) // NOLINT
    {
        const int dim = body.principal_rotation().rows();
        Eigen::VectorXd w(dim * dim);
        for (int j = 0; j < dim; ++j) {
            w.segment(dim * j, dim).array() =
                body.second_moment().diagonal()(j);
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
    const std::vector<rigid::Pose>& targets,
    const double multiplier_scaling)
{
    assert(targets.size() == bodies.num_bodies());
    assert(multiplier_scaling > 0);
    m_dim = bodies.dim();
    assert(
        x0.size()
        == (m_dim + m_dim * m_dim) * Eigen::Index(bodies.num_bodies()));

    m_x0 = x0;
    m_targets = targets;
    m_multiplier_scaling = multiplier_scaling;

    // The multipliers carry over from the previous step as the warm start,
    // unless the scene changed shape (first step or a different body count).
    if (m_lambda.size() != targets.size()) {
        m_lambda.assign(targets.size(), VectorMax3d::Zero(m_dim));
        m_lambda_A.assign(targets.size(), Eigen::VectorXd::Zero(m_dim * m_dim));
    }

    m_target_rotations.resize(targets.size());
    m_angular_weights.resize(targets.size());
    m_num_kinematic = 0;
    for (size_t i = 0; i < targets.size(); ++i) {
        // Only kinematic bodies have a meaningful target rotation; the rest
        // are left default (never read — see the KINEMATIC guards below).
        if (bodies[i].type() == rigid::RigidBody::Type::KINEMATIC) {
            m_target_rotations[i] = targets[i].rotation_matrix();
            m_angular_weights[i] = angular_weights(bodies[i]);
            ++m_num_kinematic;

            // A rotation constraint reacts with a torque, which lies in the
            // tangent space {S R̂ : S skew} of SO(dim) at the target. The rest
            // of λ_A is pushback against leaving SO(dim) that the
            // orthogonality potential provides, so we keep only the torque.
            // Carrying the rest over would push A off SO(dim) a little more
            // every step.
            const MatrixMax3d& R = m_target_rotations[i];
            const MatrixMax3d G = m_lambda_A[i].reshaped(m_dim, m_dim);
            const MatrixMax3d S = 0.5 * (G * R.transpose() - R * G.transpose());
            m_lambda_A[i] = (S * R).reshaped();
        } else {
            // A body that stopped being kinematic (e.g., converted to STATIC)
            // has no constraint, so it carries no reaction into later steps.
            reset_multipliers(i);
        }
    }

    m_kappa_q = m_params.initial_penalty;
    m_kappa_A = m_params.initial_penalty;

    m_linear_satisfied = m_num_kinematic == 0;
    m_angular_satisfied = m_num_kinematic == 0;

    m_prev_eta_q = -std::numeric_limits<double>::infinity();
    m_prev_eta_A = -std::numeric_limits<double>::infinity();
}

void AugmentedLagrangian::reset_multipliers(const size_t body_id)
{
    if (body_id >= m_lambda.size()) {
        return; // Not initialized yet, so there is nothing to carry over.
    }
    m_lambda[body_id] = VectorMax3d::Zero(m_dim);
    m_lambda_A[body_id] = Eigen::VectorXd::Zero(m_dim * m_dim);
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
    // At an energy-converged x the penalty and multiplier together supply the
    // constraint force, so folding the penalty's share into the multiplier
    // leaves λ equal to the full reaction. We do this both on the usual
    // multiplier branch and once more when a channel is satisfied, so the
    // next step starts from this step's reaction.
    //
    // Near the target, the multiplier update alone normally closes the gap
    // in one solve. It cannot when the load on the body moves between solves
    // (e.g., friction re-lagged in between), because λ then chases a moving
    // reaction and η can hover just short of satisfied indefinitely. So when
    // the gap stops shrinking, we also double the penalty, which makes the
    // body proportionally less sensitive to the change in load.
    const auto stagnated = [&](const double eta, const double prev_eta) {
        return 1 - eta > m_params.min_gap_reduction * (1 - prev_eta);
    };

    if (!m_linear_satisfied) {
        const double eta_q = linear_progress(bodies, x);
        if (eta_q >= m_params.satisfied_progress) {
            m_linear_satisfied = true;
            update_linear_multipliers(bodies, x);
        } else if (
            eta_q < m_params.stall_progress
            && m_kappa_q < m_params.max_penalty) {
            m_kappa_q *= 2;
        } else {
            update_linear_multipliers(bodies, x);
            if (stagnated(eta_q, m_prev_eta_q)
                && m_kappa_q < m_params.max_penalty) {
                m_kappa_q *= 2;
            }
        }
        m_prev_eta_q = eta_q;
        logger().debug(
            "augmented lagrangian (linear): η_q={:g} κ_q={:g} satisfied={}",
            eta_q, m_kappa_q, m_linear_satisfied);
    }

    if (!m_angular_satisfied) {
        const double eta_A = angular_progress(bodies, x); // NOLINT
        if (eta_A >= m_params.satisfied_progress) {
            m_angular_satisfied = true;
            update_angular_multipliers(bodies, x);
        } else if (
            eta_A < m_params.stall_progress
            && m_kappa_A < m_params.max_penalty) {
            m_kappa_A *= 2;
        } else {
            update_angular_multipliers(bodies, x);
            if (stagnated(eta_A, m_prev_eta_A)
                && m_kappa_A < m_params.max_penalty) {
                m_kappa_A *= 2;
            }
        }
        m_prev_eta_A = eta_A;
        logger().debug(
            "augmented lagrangian (angular): η_A={:g} κ_A={:g} satisfied={}",
            eta_A, m_kappa_A, m_angular_satisfied);
    }
}

void AugmentedLagrangian::update_linear_multipliers(
    const rigid::RigidBodies& bodies, Eigen::ConstRef<Eigen::VectorXd> x)
{
    const int ndof = m_dim + m_dim * m_dim;
    for (size_t i = 0; i < bodies.num_bodies(); ++i) {
        if (bodies[i].type() != rigid::RigidBody::Type::KINEMATIC) {
            continue;
        }
        m_lambda[i] -= m_kappa_q * bodies[i].mass()
            * (x.segment(ndof * i, m_dim) - m_targets[i].position)
            / m_multiplier_scaling;
    }
}

void AugmentedLagrangian::update_angular_multipliers(
    const rigid::RigidBodies& bodies, Eigen::ConstRef<Eigen::VectorXd> x)
{
    const int ndof = m_dim + m_dim * m_dim;
    for (size_t i = 0; i < bodies.num_bodies(); ++i) {
        if (bodies[i].type() != rigid::RigidBody::Type::KINEMATIC) {
            continue;
        }
        const Eigen::VectorXd da = x.segment(ndof * i + m_dim, m_dim * m_dim)
            - Eigen::VectorXd(m_target_rotations[i].reshaped());
        m_lambda_A[i] -= m_kappa_A * m_angular_weights[i].cwiseProduct(da)
            / m_multiplier_scaling;
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
            - m_multiplier_scaling * m_lambda[i].dot(dq);

        const Eigen::VectorXd da = x.segment(ndof * i + m_dim, m_dim * m_dim)
            - Eigen::VectorXd(m_target_rotations[i].reshaped());
        const Eigen::VectorXd& w = m_angular_weights[i];
        energy += 0.5 * m_kappa_A * da.dot(w.cwiseProduct(da))
            - m_multiplier_scaling * m_lambda_A[i].dot(da);
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
            - m_multiplier_scaling * m_lambda[i];

        const Eigen::VectorXd da = x.segment(ndof * i + m_dim, m_dim * m_dim)
            - Eigen::VectorXd(m_target_rotations[i].reshaped());
        const Eigen::VectorXd& w = m_angular_weights[i];
        grad.segment(ndof * i + m_dim, m_dim * m_dim) =
            m_kappa_A * w.cwiseProduct(da)
            - m_multiplier_scaling * m_lambda_A[i];
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
        const Eigen::VectorXd& w = m_angular_weights[i];
        for (int k = 0; k < m_dim * m_dim; ++k) {
            triplets.emplace_back(
                ndof * i + m_dim + k, ndof * i + m_dim + k, m_kappa_A * w(k));
        }
    }
    hess.setFromTriplets(triplets.begin(), triplets.end());
    return hess;
}

} // namespace ipc::affine
