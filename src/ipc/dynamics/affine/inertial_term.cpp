#include "inertial_term.hpp"

#include <ipc/dynamics/time_integration/time_integrator.hpp>

namespace ipc::affine {

InertialTerm::InertialTerm(
    const rigid::RigidBodies& bodies,
    const std::shared_ptr<const dynamics::ImplicitTimeIntegrator>&
        _time_integrator)
    : time_integrator(_time_integrator)
{
    build_mass_matrix(bodies);
}

void InertialTerm::build_mass_matrix(const rigid::RigidBodies& bodies)
{
    // The exact ABD mass matrix: because the rest positions are COM-centered
    // in the principal inertia frame, ∫ρ x̄ dV = 0 and ∫ρ x̄x̄ᵀ dV = J
    // (diagonal), so with q = [p; vec(A) column-major] the inertial energy
    // ½m‖p − p̂‖² + ½tr((A − Â) J (A − Â)ᵀ) has mass matrix
    // M = blkdiag(m I, J ⊗ I).
    // Non-DYNAMIC bodies contribute no inertia (zero blocks): their motion is
    // fully prescribed (pinned DOFs and/or the augmented Lagrangian).
    const int dim = bodies.dim();
    const int ndof = dim + dim * dim;
    std::vector<Eigen::Triplet<double>> triplets;
    for (int i = 0; i < bodies.num_bodies(); i++) {
        if (!bodies[i].is_dynamic()) {
            continue;
        }
        const double mass = bodies[i].mass();
        const auto& J = bodies[i].J().diagonal();
        assert(J.size() == dim);
        for (int k = 0; k < dim; k++) {
            triplets.emplace_back(i * ndof + k, i * ndof + k, mass);
            for (int j = 0; j < dim; j++) {
                const int dof = i * ndof + dim + k + dim * j; // A(k, j)
                triplets.emplace_back(dof, dof, J(j));
            }
        }
    }
    m_mass.resize(ndof * bodies.num_bodies(), ndof * bodies.num_bodies());
    m_mass.setFromTriplets(triplets.begin(), triplets.end());
}

void InertialTerm::update(const rigid::RigidBodies& bodies)
{
    // Update the predicted positions from the current time integrator state
    m_x_hat = time_integrator->predicted_positions();
    // Body types can change at runtime (convert_to_static), so refresh the
    // zero blocks of non-DYNAMIC bodies.
    build_mass_matrix(bodies);
}

// ---- Cumulative functions ---------------------------------------------------

double InertialTerm::operator()(
    const rigid::RigidBodies& bodies, Eigen::ConstRef<Eigen::VectorXd> x) const
{
    assert(x.size() == m_x_hat.size());
    Eigen::VectorXd dx = x - m_x_hat;
    return 0.5 * (dx.transpose() * m_mass * dx)(0, 0);
}

Eigen::VectorXd InertialTerm::gradient(
    const rigid::RigidBodies& bodies, Eigen::ConstRef<Eigen::VectorXd> x) const
{
    assert(x.size() == m_x_hat.size());
    return m_mass * (x - m_x_hat);
}

} // namespace ipc::affine