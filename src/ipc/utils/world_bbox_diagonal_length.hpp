#pragma once

#include <Eigen/Core>

namespace ipc {

inline double world_bbox_diagonal_length(const Eigen::MatrixXd& V)
{
    return (V.colwise().maxCoeff() - V.colwise().minCoeff()).norm();
}

} // namespace ipc
