#pragma once

#include <ipc/utils/eigen_ext.hpp>

#include <array>
#include <limits>

namespace ipc {

class CollisionStencil {
public:
    virtual ~CollisionStencil() = default;

    /// @brief Get the number of vertices in the contact stencil.
    virtual int num_vertices() const = 0;

    /// @brief Get the vertex IDs of the contact stencil.
    /// @param edges Collision mesh edges
    /// @param faces Collision mesh faces
    /// @return The vertex IDs of the contact stencil. Size is always 4, but elements i > num_vertices() are -1.
    virtual std::array<long, 4> vertex_ids(
        const Eigen::MatrixXi& edges, const Eigen::MatrixXi& faces) const = 0;

    /// @brief Get the vertex positions of the contact stencil.
    /// @param vertices Vertex positions
    /// @param edges Collision mesh edges
    /// @param faces Collision mesh faces
    /// @return The vertex positions of the contact stencil. Size is always 4, but elements i > num_vertices() are NaN.
    std::array<VectorMax3d, 4> vertices(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces) const;
};

} // namespace ipc