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

Eigen::SparseMatrix<double> RigidPotential::hessian(
    const RigidBodies& bodies,
    Eigen::ConstRef<Eigen::VectorXd> x,
    const PSDProjectionMethod project_hessian_to_psd) const
{
    const int ndof = x.size() / bodies.num_bodies();

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(bodies.num_bodies() * ndof * ndof);
    for (size_t i = 0; i < bodies.num_bodies(); ++i) {
        const MatrixMax6d block = hessian(
            i, bodies[i], x.segment(i * ndof, ndof), project_hessian_to_psd);
        assert(block.allFinite());
        for (int c = 0; c < ndof; ++c) {
            for (int r = 0; r < ndof; ++r) {
                if (block(r, c) != 0) {
                    triplets.emplace_back(
                        i * ndof + r, i * ndof + c, block(r, c));
                }
            }
        }
    }

    Eigen::SparseMatrix<double> hess(x.size(), x.size());
    hess.setFromTriplets(triplets.begin(), triplets.end());
    return hess;
}

} // namespace ipc::rigid