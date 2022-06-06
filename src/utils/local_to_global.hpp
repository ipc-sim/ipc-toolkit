#pragma once

#include <vector>

#include <Eigen/Core>
#include <Eigen/Sparse>

namespace ipc {

// IDs can be either std::array<long, N> or std::vector<long>
template <typename DerivedLocalGrad, typename IDs, typename DerivedGrad>
void local_gradient_to_global_gradient(
    const Eigen::MatrixBase<DerivedLocalGrad>& local_grad,
    const IDs& ids,
    int dim,
    Eigen::PlainObjectBase<DerivedGrad>& grad)
{
    assert(local_grad.size() % dim == 0);
    const int n_verts = local_grad.size() / dim;
    assert(ids.size() >= n_verts); // Can be extra ids
    for (int i = 0; i < n_verts; i++) {
        grad.segment(dim * ids[i], dim) += local_grad.segment(dim * i, dim);
    }
}

// IDs can be either std::array<long, N> or std::vector<long>
template <typename Derived, typename IDs>
void local_hessian_to_global_triplets(
    const Eigen::MatrixBase<Derived>& local_hessian,
    const IDs& ids,
    int dim,
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
