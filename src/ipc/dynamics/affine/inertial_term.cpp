#include "inertial_term.hpp"

#include <ipc/dynamics/rigid/time_integrator.hpp>

namespace ipc::affine {

InertialTerm::InertialTerm(
    const rigid::RigidBodies& bodies,
    const std::shared_ptr<const rigid::ImplicitEuler>& _time_integrator)
    : time_integrator(_time_integrator)
{
    std::vector<Eigen::Triplet<double>> triplets;
    for (int i = 0; i < bodies.num_bodies(); i++) {
        const Eigen::SparseMatrix<double> J =
            rigid::AffinePose::J(bodies.body_rest_positions(i));
        const MatrixMax12d Mi = bodies[i].density() * (J.transpose() * J);
        for (int j = 0; j < Mi.rows(); j++) {
            for (int k = 0; k < Mi.cols(); k++) {
                triplets.emplace_back(i * 12 + j, i * 12 + k, Mi(j, k));
            }
        }
    }
    m_mass.resize(12 * bodies.num_bodies(), 12 * bodies.num_bodies());
    m_mass.setFromTriplets(triplets.begin(), triplets.end());
}

void InertialTerm::update(const rigid::RigidBodies& bodies)
{
    // Update the predicted poses based on the current time integrator state
    m_x_hat = time_integrator->x_hat();
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

Eigen::MatrixXd InertialTerm::hessian(
    const rigid::RigidBodies& bodies,
    Eigen::ConstRef<Eigen::VectorXd> x,
    const PSDProjectionMethod project_hessian_to_psd) const
{
    return m_mass;
}

} // namespace ipc::affine