#pragma once

#include <vector>

#include <Eigen/Core>
#include <Eigen/Sparse>

namespace ipc {

void local_hessian_to_global_triplets(
    const Eigen::MatrixXd& local_hessian,
    const std::vector<long>& ids,
    int dim,
    std::vector<Eigen::Triplet<double>>& triplets);

} // namespace ipc
