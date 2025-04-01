#pragma once

#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

/// @brief Compute the diagonal length of the world bounding box.
/// @param vertices Vertex positions
/// @return The diagonal length of the world bounding box.
inline double
world_bbox_diagonal_length(Eigen::ConstRef<Eigen::MatrixXd> vertices)
{
    return (vertices.colwise().maxCoeff() - vertices.colwise().minCoeff())
        .norm();
}

} // namespace ipc
