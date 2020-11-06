#pragma once

#include <vector>

#include <Eigen/Core>
#include <Eigen/Sparse>

namespace ipc {

template <typename Derived>
void local_hessian_to_global_triplets(
    const Eigen::MatrixBase<Derived>& local_hessian,
    const std::vector<long>& ids,
    int dim,
    std::vector<Eigen::Triplet<double>>& triplets)
{
    for (int i = 0; i < ids.size(); i++) {
        for (int j = 0; j < ids.size(); j++) {
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
