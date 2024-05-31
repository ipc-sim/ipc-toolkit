#include "orthogonality_potential.hpp"

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>

namespace ipc::affine {

double
OrthogonalityPotential::operator()(const std::vector<AffineBody>& bodies) const
{
    if (bodies.empty()) {
        return 0;
    }

    tbb::enumerable_thread_specific<double> storage(0);

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), bodies.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            double& local_potential = storage.local();
            for (size_t i = r.begin(); i < r.end(); ++i) {
                // Quadrature weight is premultiplied by local potential
                local_potential += (*this)(bodies[i]);
            }
        });

    return storage.combine([](double a, double b) { return a + b; });
}

Eigen::VectorXd
OrthogonalityPotential::gradient(const std::vector<AffineBody>& bodies) const
{
    if (bodies.empty()) {
        return Eigen::VectorXd();
    }

    const int dim = bodies[0].A.rows();
    const int ndof_per_body = dim * dim + dim;
    Eigen::VectorXd grad(ndof_per_body * bodies.size());

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), bodies.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i < r.end(); ++i) {
                grad.segment(ndof_per_body * i, ndof_per_body) =
                    this->gradient(bodies[i]);
            }
        });

    return grad;
}

Eigen::SparseMatrix<double> OrthogonalityPotential::hessian(
    const std::vector<AffineBody>& bodies,
    const bool project_hessian_to_psd) const
{
    if (bodies.empty()) {
        return Eigen::SparseMatrix<double>();
    }

    const int dim = bodies[0].A.rows();
    const int ndof_per_body = dim * dim + dim;
    const int ndof = ndof_per_body * bodies.size();

    tbb::enumerable_thread_specific<std::vector<Eigen::Triplet<double>>>
        storage;

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), bodies.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            auto& hess_triplets = storage.local();

            for (size_t i = r.begin(); i < r.end(); ++i) {
                const MatrixMax9d local_hess = this->hessian(bodies[i]);
                for (size_t hi = 0; hi < local_hess.rows(); ++hi) {
                    for (size_t hj = 0; hj < local_hess.rows(); ++hj) {
                        hess_triplets.emplace_back(
                            ndof_per_body * i + hi, ndof_per_body * i + hj,
                            local_hess(hi, hj));
                    }
                }
            }
        });

    Eigen::SparseMatrix<double> hess(ndof, ndof);
    for (const auto& local_hess_triplets : storage) {
        Eigen::SparseMatrix<double> local_hess(ndof, ndof);
        local_hess.setFromTriplets(
            local_hess_triplets.begin(), local_hess_triplets.end());
        hess += local_hess;
    }
    return hess;
}

// -- Single body methods ---------------------------------------------

double OrthogonalityPotential::operator()(const AffineBody& body) const
{
    const auto& [A, p, volume] = body;
    const auto I = MatrixMax3d::Identity(A.rows(), A.cols());
    return stiffness * volume * (A * A.transpose() - I).squaredNorm();

    // double r = 0;
    // for (int i = 0; i < A.cols(); i++) {
    //     for (int j = 0; j < A.cols(); j++) {
    //         const double dot = A.col(i).dot(A.col(j)) - int(i == j);
    //         r += dot * dot;
    //     }
    // }
    // return stiffness * volume * r;
}

VectorMax12d OrthogonalityPotential::gradient(const AffineBody& body) const
{
    const auto& [A, p, volume] = body;
    VectorMax12d grad = VectorMax12d::Zero(A.size() + p.size());

    const auto I = MatrixMax3d::Identity(A.rows(), A.cols());
    const MatrixMax3d G = stiffness * volume * 4 * (A * A.transpose() - I) * A;
    grad.tail(A.size()) = G.reshaped();

    // for (int i = 0; i < A.cols(); i++) {
    //     auto grad_ai = grad.segment(i * A.rows() + p.size(), A.rows());
    //     for (int j = 0; j < A.cols(); j++) {
    //         grad_ai += (A.col(i).dot(A.col(j)) - int(i == j)) * A.col(j);
    //     }
    //     grad_ai *= 4 * stiffness * volume;
    // }

    return grad;
}

MatrixMax12d OrthogonalityPotential::hessian(const AffineBody& body) const
{
    const auto& [A, p, volume] = body;
    const int dim = A.rows();
    assert(p.size() == dim);
    const int ndof = A.size() + p.size();
    const auto I = MatrixMax3d::Identity(A.rows(), A.cols());

    MatrixMax12d hess = MatrixMax12d::Zero(ndof, ndof);
    // NOTE: top left block is with respect to p (translation) and is zero

    for (int i = 0; i < A.cols(); i++) {
        for (int j = 0; j < A.cols(); j++) {
            auto hess_aij =
                hess.block(i * dim + p.size(), j * dim + p.size(), dim, dim);

            hess_aij += (int(i == j) + 1) * A.col(j) * A.col(i).transpose()
                + (A.col(i).dot(A.col(j)) - int(i == j)) * I;

            if (i == j) {
                for (int k = 0; k < A.cols(); k++) {
                    if (i != k) {
                        hess_aij += A.col(k) * A.col(k).transpose();
                    }
                }
            }

            hess_aij *= 4 * stiffness * volume;
        }
    }

    return hess;
}

} // namespace ipc::affine