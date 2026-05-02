#include "rigid_potential.hpp"

namespace ipc::rigid {

// ---- Cumulative functions ---------------------------------------------------

double RigidPotential::operator()(
    const RigidBodies& bodies, Eigen::ConstRef<Eigen::VectorXd> x) const
{
    const int ndof = x.size() / bodies.num_bodies();

    double energy = 0.0;
    for (size_t i = 0; i < bodies.num_bodies(); ++i) {
        energy += operator()(i, bodies[i], x.segment(i * ndof, ndof));
    }
    return energy;
}

Eigen::VectorXd RigidPotential::gradient(
    const RigidBodies& bodies, Eigen::ConstRef<Eigen::VectorXd> x) const
{
    const int ndof = x.size() / bodies.num_bodies();

    Eigen::VectorXd grad(x.size());
    for (size_t i = 0; i < bodies.num_bodies(); ++i) {
        grad.segment(i * ndof, ndof) =
            gradient(i, bodies[i], x.segment(i * ndof, ndof));
    }
    return grad;
}

Eigen::MatrixXd RigidPotential::hessian(
    const RigidBodies& bodies,
    Eigen::ConstRef<Eigen::VectorXd> x,
    const PSDProjectionMethod project_hessian_to_psd) const
{
    const int ndof = x.size() / bodies.num_bodies();

    Eigen::MatrixXd hess = Eigen::MatrixXd::Zero(x.size(), x.size());
    for (size_t i = 0; i < bodies.num_bodies(); ++i) {
        hess.block(i * ndof, i * ndof, ndof, ndof) = hessian(
            i, bodies[i], x.segment(i * ndof, ndof), project_hessian_to_psd);
    }
    assert(hess.allFinite());
    return hess;
}

} // namespace ipc::rigid