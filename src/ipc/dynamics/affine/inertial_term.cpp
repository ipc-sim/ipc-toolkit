#include "inertial_term.hpp"

#include <ipc/dynamics/time_integration/time_integrator.hpp>

namespace ipc::affine {

InertialTerm::InertialTerm(
    const rigid::RigidBodies& bodies,
    const std::shared_ptr<const dynamics::ImplicitTimeIntegrator>&
        _time_integrator)
    : time_integrator(_time_integrator)
{
    assert(bodies.dim() == 3); // TODO: Support 2D affine bodies

    // The exact ABD mass matrix: because the rest positions are COM-centered
    // in the principal inertia frame, ∫ρ x̄ dV = 0 and ∫ρ x̄x̄ᵀ dV = J
    // (diagonal), so with q = [p; vec(A) column-major] the inertial energy
    // ½m‖p − p̂‖² + ½tr((A − Â) J (A − Â)ᵀ) has mass matrix
    // M = blkdiag(m I₃, J ⊗ I₃).
    std::vector<Eigen::Triplet<double>> triplets;
    for (int i = 0; i < bodies.num_bodies(); i++) {
        const double mass = bodies[i].mass();
        const auto& J = bodies[i].J().diagonal();
        for (int k = 0; k < 3; k++) {
            triplets.emplace_back(i * 12 + k, i * 12 + k, mass);
            for (int j = 0; j < 3; j++) {
                const int dof = i * 12 + 3 + k + 3 * j; // A(k, j)
                triplets.emplace_back(dof, dof, J(j));
            }
        }
    }
    m_mass.resize(12 * bodies.num_bodies(), 12 * bodies.num_bodies());
    m_mass.setFromTriplets(triplets.begin(), triplets.end());
}

void InertialTerm::update(const rigid::RigidBodies& bodies)
{
    // Update the predicted positions from the current time integrator state
    m_x_hat = time_integrator->predicted_positions();
}

// ---- Cumulative functions ---------------------------------------------------

double InertialTerm::operator()(
    const rigid::RigidBodies& bodies, Eigen::ConstRef<Eigen::VectorXd> x) const
{
    Eigen::VectorXd dx = x - m_x_hat;
    return 0.5 * (dx.transpose() * m_mass * dx)(0, 0);
}

Eigen::VectorXd InertialTerm::gradient(
    const rigid::RigidBodies& bodies, Eigen::ConstRef<Eigen::VectorXd> x) const
{
    return m_mass * (x - m_x_hat);
}

} // namespace ipc::affine