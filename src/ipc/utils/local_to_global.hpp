#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <vector>

namespace ipc {

template <typename IDContainer>
void local_gradient_to_global_gradient(
    Eigen::ConstRef<Eigen::VectorXd> local_grad,
    const IDContainer& ids,
    const int dim,
    Eigen::Ref<Eigen::VectorXd> grad)
{
    assert(local_grad.size() % dim == 0);
    const int n_verts = local_grad.size() / dim;
    assert(ids.size() >= n_verts); // Can be extra ids
    for (int i = 0; i < n_verts; i++) {
        grad.segment(dim * ids[i], dim) += local_grad.segment(dim * i, dim);
    }
}

template <typename IDContainer>
void local_gradient_to_global_gradient(
    Eigen::ConstRef<Eigen::VectorXd> local_grad,
    const IDContainer& ids,
    const int dim,
    Eigen::SparseVector<double>& grad)
{
    assert(local_grad.size() % dim == 0);
    const int n_verts = local_grad.size() / dim;
    assert(ids.size() >= n_verts); // Can be extra ids
    for (int i = 0; i < n_verts; i++) {
        for (int d = 0; d < dim; d++) {
            grad.coeffRef(dim * ids[i] + d) += local_grad(dim * i + d);
        }
    }
}

template <typename IDContainer>
void local_hessian_to_global_triplets(
    Eigen::ConstRef<Eigen::MatrixXd> local_hessian,
    const IDContainer& ids,
    const int dim,
    std::vector<Eigen::Triplet<double>>& triplets)
{
    assert(local_hessian.rows() == local_hessian.cols());
    assert(local_hessian.rows() % dim == 0);
    const int n_verts = local_hessian.rows() / dim;
    assert(ids.size() >= n_verts); // Can be extra ids
    for (int i = 0; i < n_verts; i++) {
        for (int j = 0; j < n_verts; j++) {
            for (int k = 0; k < dim; k++) {
                for (int l = 0; l < dim; l++) {
                    triplets.emplace_back(
                        dim * ids[i] + k, dim * ids[j] + l,
                        local_hessian(dim * i + k, dim * j + l));
                }
            }
        }
    }
}

} // namespace ipc
