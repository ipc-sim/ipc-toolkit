#pragma once

#include <ipc/dynamics/rigid/rigid_body.hpp>

#include <memory>

namespace ipc::rigid {

class RigidBodies : public CollisionMesh {
    RigidBodies() = delete;

public:
    RigidBodies(
        Eigen::ConstRef<Eigen::MatrixXd> rest_positions,
        Eigen::ConstRef<Eigen::MatrixXi> edges,
        Eigen::ConstRef<Eigen::MatrixXi> faces,
        const std::vector<index_t>& body_vertex_starts,
        const std::vector<index_t>& body_edge_starts,
        const std::vector<index_t>& body_face_starts,
        const std::vector<double>& densities,
        std::vector<Pose>& initial_poses);

    static std::shared_ptr<RigidBodies> build_from_meshes(
        const std::vector<Eigen::MatrixXd>& rest_positions,
        const std::vector<Eigen::MatrixXi>& edges,
        const std::vector<Eigen::MatrixXi>& faces,
        const std::vector<double>& densities,
        std::vector<Pose>& initial_poses);

    Eigen::MatrixXd vertices(const std::vector<Pose>& poses) const
    {
        assert(poses.size() == num_bodies());
        Eigen::MatrixXd V(num_vertices(), dim());
        for (size_t i = 0; i < num_bodies(); ++i) {
            const index_t start = body_vertex_starts[i];
            const index_t end = body_vertex_starts[i + 1];
            V.middleRows(start, end - start) = poses[i].transform_vertices(
                rest_positions().middleRows(
                    body_vertex_starts[i], end - start));
        }
        return V;
    }

    Eigen::MatrixXd vertices(
        Eigen::ConstRef<Eigen::MatrixXd> full_positions_or_poses) const override
    {
        if (full_positions_or_poses.cols() == dim()) {
            assert(full_positions_or_poses.rows() == full_num_vertices());
            return CollisionMesh::vertices(full_positions_or_poses);
        } else {
            assert(full_positions_or_poses.cols() == 1);
            return vertices(Pose::to_poses(full_positions_or_poses, dim()));
        }
    }

    /// @brief Get the rigid body at index i.
    /// @param i Index of the rigid body.
    /// @return Reference to the rigid body at index i.
    RigidBody& operator[](size_t i) { return bodies[i]; }

    /// @brief Get the rigid body at index i.
    /// @param i Index of the rigid body.
    /// @return Const reference to the rigid body at index i.
    const RigidBody& operator[](size_t i) const { return bodies[i]; }

    /// @brief Get the number of rigid bodies in the system.
    /// @return Number of rigid bodies.
    size_t num_bodies() const { return bodies.size(); }

    /// @brief Get the number of vertices in the i-th rigid body mesh.
    /// @param i Index of the rigid body mesh.
    /// @return Number of vertices in the i-th rigid body mesh.
    size_t body_num_vertices(size_t i) const
    {
        return body_vertex_starts[i + 1] - body_vertex_starts[i];
    }

    /// @brief Get the number of edges in the i-th rigid body mesh.
    /// @param i Index of the rigid body mesh.
    /// @return Number of edges in the i-th rigid body mesh.
    size_t body_num_edges(size_t i) const
    {
        return body_edge_starts[i + 1] - body_edge_starts[i];
    }

    /// @brief Get the number of faces in the i-th rigid body mesh.
    /// @param i Index of the rigid body mesh.
    /// @return Number of faces in the i-th rigid body mesh.
    size_t body_num_faces(size_t i) const
    {
        return body_face_starts[i + 1] - body_face_starts[i];
    }

    /// @brief Get the vertices of the i-th rigid body mesh.
    /// @param i Index of the rigid body mesh.
    /// @return Vertices of the i-th rigid body mesh.
    auto body_rest_positions(size_t i) const
    {
        return rest_positions().middleRows(
            body_vertex_starts[i], body_num_vertices(i));
    }

    /// @brief Get the vertices of the i-th rigid body mesh.
    /// @param i Index of the rigid body mesh.
    /// @return Vertices of the i-th rigid body mesh.
    Eigen::MatrixXd body_vertices(size_t i, const Pose& pose) const
    {
        return pose.transform_vertices(body_rest_positions(i));
    }

    /// @brief Get the edges of the i-th rigid body mesh.
    /// @note Returns indices in the local body vertex indexing.
    /// @param i Index of the rigid body mesh.
    /// @return Edges of the i-th rigid body mesh.
    Eigen::MatrixXi body_edges(size_t i) const
    {
        return edges()
                   .middleRows(body_edge_starts[i], body_num_edges(i))
                   .array()
            - body_vertex_starts[i];
    }

    /// @brief Get the faces of the i-th rigid body mesh.
    /// @note Returns indices in the local body vertex indexing.
    /// @param i Index of the rigid body mesh.
    /// @return Faces of the i-th rigid body mesh.
    Eigen::MatrixXi body_faces(size_t i) const
    {
        return faces()
                   .middleRows(body_face_starts[i], body_num_faces(i))
                   .array()
            - body_vertex_starts[i];
    }

    /// @brief Map a vector quantity on the collision mesh to the rigid degrees of freedom.
    /// This is useful for mapping gradients from the collision mesh to the
    /// rigid degrees of freedom (i.e., applies the chain-rule).
    /// @param x Vector quantity on the collision mesh with size equal to ndof().
    /// @return Vector quantity on the full mesh with size equal to full_ndof().
    Eigen::VectorXd to_rigid_dof(
        const std::vector<Pose>& poses,
        Eigen::ConstRef<Eigen::VectorXd> g) const;

    /// @brief Map a matrix quantity on the collision mesh to the full mesh.
    /// This is useful for mapping Hessians from the collision mesh to the full
    /// mesh (i.e., applies the chain-rule).
    /// @param X Matrix quantity on the collision mesh with size equal to ndof() × ndof().
    /// @return Matrix quantity on the full mesh with size equal to full_ndof() × full_ndof().
    Eigen::MatrixXd to_rigid_dof(
        const std::vector<Pose>& poses,
        Eigen::ConstRef<Eigen::VectorXd> g,
        const Eigen::SparseMatrix<double>& H) const;

private:
    /// @brief Get the body index for a given vertex index.
    index_t vertex_to_body(index_t vi) const
    {
        assert(vi >= 0 && vi < num_vertices());

        // Find the first element GREATER than vi
        auto it = std::upper_bound(
            body_vertex_starts.begin(), body_vertex_starts.end(), vi);

        // The body index is (iterator distance - 1)
        // std::distance returns the number of elements from the start to 'it'
        return static_cast<index_t>(
                   std::distance(body_vertex_starts.begin(), it))
            - 1;
    }

    /// @brief Rigid bodies in the system
    std::vector<RigidBody> bodies;

    // Start indices of vertices for each body
    std::vector<index_t> body_vertex_starts;
    std::vector<index_t> body_edge_starts;
    std::vector<index_t> body_face_starts;
};

} // namespace ipc::rigid