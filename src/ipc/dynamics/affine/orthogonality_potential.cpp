#include "orthogonality_potential.hpp"

#include <oneapi/tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>

namespace ipc::affine {

double OrthogonalityPotential::operator()(
    const rigid::RigidBodies& bodies, Eigen::ConstRef<Eigen::VectorXd> x) const
{
    if (bodies.num_bodies() == 0) {
        return 0;
    }

    const int ndof = x.size() / bodies.num_bodies();

    return tbb::parallel_reduce(
        tbb::blocked_range<size_t>(size_t(0), bodies.num_bodies()), 0.0,
        [&](const tbb::blocked_range<size_t>& r, double local_sum) {
            for (size_t i = r.begin(); i < r.end(); ++i) {
                local_sum += operator()(bodies[i], x.segment(i * ndof, ndof));
            }
            return local_sum;
        },
        std::plus<double>());
}

Eigen::VectorXd OrthogonalityPotential::gradient(
    const rigid::RigidBodies& bodies, Eigen::ConstRef<Eigen::VectorXd> x) const
{
    if (bodies.num_bodies() == 0) {
        return Eigen::VectorXd();
    }

    const int ndof = x.size() / bodies.num_bodies();

    Eigen::VectorXd grad(x.size());
    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), bodies.num_bodies()),
        [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i < r.end(); ++i) {
                grad.segment(ndof * i, ndof) =
                    gradient(bodies[i], x.segment(i * ndof, ndof));
            }
        });
    return grad;
}

Eigen::SparseMatrix<double> OrthogonalityPotential::hessian(
    const rigid::RigidBodies& bodies,
    Eigen::ConstRef<Eigen::VectorXd> x,
    const PSDProjectionMethod project_hessian_to_psd) const
{
    if (bodies.num_bodies() == 0) {
        return Eigen::SparseMatrix<double>();
    }

    const int ndof = x.size() / bodies.num_bodies();

    tbb::enumerable_thread_specific<std::vector<Eigen::Triplet<double>>>
        storage;

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), bodies.num_bodies()),
        [&](const tbb::blocked_range<size_t>& r) {
            auto& hess_triplets = storage.local();

            for (size_t i = r.begin(); i < r.end(); ++i) {
                const MatrixMax12d local_hess = hessian(
                    bodies[i], x.segment(i * ndof, ndof),
                    project_hessian_to_psd);
                for (size_t hi = 0; hi < local_hess.rows(); ++hi) {
                    for (size_t hj = 0; hj < local_hess.rows(); ++hj) {
                        hess_triplets.emplace_back(
                            ndof * i + hi, ndof * i + hj, local_hess(hi, hj));
                    }
                }
            }
        });

    Eigen::SparseMatrix<double> hess(x.size(), x.size());
    for (const auto& local_hess_triplets : storage) {
        Eigen::SparseMatrix<double> local_hess(x.size(), x.size());
        local_hess.setFromTriplets(
            local_hess_triplets.begin(), local_hess_triplets.end());
        hess += local_hess;
    }
    return hess;
}

// -- Single body methods ---------------------------------------------

double OrthogonalityPotential::operator()(
    const rigid::RigidBody& body, Eigen::ConstRef<VectorMax12d> x) const
{
    const int dim = body.R0().rows();
    const auto& A = x.tail(dim * dim).reshaped(dim, dim);
    const auto I = MatrixMax3d::Identity(dim, dim);
    return stiffness * body.volume() * (A * A.transpose() - I).squaredNorm();
}

VectorMax12d OrthogonalityPotential::gradient(
    const rigid::RigidBody& body, Eigen::ConstRef<VectorMax12d> x) const
{
    VectorMax12d grad = VectorMax12d::Zero(x.size());

    const int dim = body.R0().rows();
    const auto& A = x.tail(dim * dim).reshaped(dim, dim);
    const auto I = MatrixMax3d::Identity(dim, dim);

    const MatrixMax3d G =
        stiffness * body.volume() * 4 * (A * A.transpose() - I) * A;
    grad.tail(A.size()) = G.reshaped();

    return grad;
}

MatrixMax12d OrthogonalityPotential::hessian(
    const rigid::RigidBody& body,
    Eigen::ConstRef<VectorMax12d> x,
    const PSDProjectionMethod project_hessian_to_psd) const
{
    const int dim = body.R0().rows();
    const auto& A = x.tail(dim * dim).reshaped(dim, dim);
    const auto I = MatrixMax3d::Identity(dim, dim);

    MatrixMax12d hess = MatrixMax12d::Zero(x.size(), x.size());
    // NOTE: top left block is with respect to p (translation) and is zero

    for (int i = 0; i < A.cols(); i++) {
        for (int j = 0; j < A.cols(); j++) {
            auto hess_aij = hess.block(i * dim + dim, j * dim + dim, dim, dim);

            // Extra outer product term if i == j
            hess_aij += (int(i == j) + 1) * A.col(j) * A.col(i).transpose()
                + (A.col(i).dot(A.col(j)) - int(i == j)) * I;

            if (i == j) {
                for (int k = 0; k < A.cols(); k++) {
                    if (i != k) { // Already added the i == k term above
                        hess_aij += A.col(k) * A.col(k).transpose();
                    }
                }
            }

            hess_aij *= 4 * stiffness * body.volume();
        }
    }

    return project_to_psd(hess, project_hessian_to_psd);
}

} // namespace ipc::affine