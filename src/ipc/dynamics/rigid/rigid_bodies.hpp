#pragma once

#include <ipc/dynamics/rigid/rigid_body.hpp>

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

    static RigidBodies build_from_meshes(
        const std::vector<Eigen::MatrixXd>& rest_positions,
        const std::vector<Eigen::MatrixXi>& edges,
        const std::vector<Eigen::MatrixXi>& faces,
        const std::vector<double>& densities,
        std::vector<Pose>& initial_poses);

    Eigen::MatrixXd vertices(const std::vector<Pose>& poses) const
    {
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

    /// @brief Get the rigid body at index i.
    /// @param i Index of the rigid body.
    /// @return Reference to the rigid body at index i.
    const RigidBody& operator[](size_t i) const { return bodies[i]; }

    /// @brief Get the number of rigid bodies in the system.
    /// @return Number of rigid bodies.
    size_t num_bodies() const { return bodies.size(); }

    /// @brief Get the number of vertices in the i-th rigid body mesh.
    /// @param i Index of the rigid body mesh.
    /// @return Number of vertices in the i-th rigid body mesh.
    size_t num_body_vertices(size_t i) const
    {
        return body_vertex_starts[i + 1] - body_vertex_starts[i];
    }

    /// @brief Get the number of edges in the i-th rigid body mesh.
    /// @param i Index of the rigid body mesh.
    /// @return Number of edges in the i-th rigid body mesh.
    size_t num_body_edges(size_t i) const
    {
        return body_edge_starts[i + 1] - body_edge_starts[i];
    }

    /// @brief Get the number of faces in the i-th rigid body mesh.
    /// @param i Index of the rigid body mesh.
    /// @return Number of faces in the i-th rigid body mesh.
    size_t num_body_faces(size_t i) const
    {
        return body_face_starts[i + 1] - body_face_starts[i];
    }

    /// @brief Get the vertices of the i-th rigid body mesh.
    /// @param i Index of the rigid body mesh.
    /// @return Vertices of the i-th rigid body mesh.
    auto body_vertices(size_t i) const
    {
        return rest_positions().middleRows(
            body_vertex_starts[i], num_body_vertices(i));
    }

    /// @brief Get the edges of the i-th rigid body mesh.
    /// @param i Index of the rigid body mesh.
    /// @return Edges of the i-th rigid body mesh.
    auto body_edges(size_t i) const
    {
        return edges().middleRows(body_edge_starts[i], num_body_edges(i));
    }

    /// @brief Get the faces of the i-th rigid body mesh.
    /// @param i Index of the rigid body mesh.
    /// @return Faces of the i-th rigid body mesh.
    auto body_faces(size_t i) const
    {
        return faces().middleRows(body_face_starts[i], num_body_faces(i));
    }

private:
    std::vector<RigidBody> bodies;

    // Start indices of vertices for each body
    std::vector<index_t> body_vertex_starts;
    std::vector<index_t> body_edge_starts;
    std::vector<index_t> body_face_starts;
};

} // namespace ipc::rigid