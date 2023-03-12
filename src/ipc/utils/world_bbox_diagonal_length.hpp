#pragma once

#include <Eigen/Core>

namespace ipc {

/// @brief Compute the diagonal length of the world bounding box.
/// @param V Positions of the vertices.
/// @return The diagonal length of the world bounding box.
inline double world_bbox_diagonal_length(const Eigen::MatrixXd& V)
{
    return (V.colwise().maxCoeff() - V.colwise().minCoeff()).norm();
}

} // namespace ipc
