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

    size_t num_bodies() const { return bodies.size(); }

    const RigidBody& operator[](size_t i) const { return bodies[i]; }

private:
    std::vector<RigidBody> bodies;

    // Start indices of vertices for each body
    std::vector<index_t> body_vertex_starts;
    std::vector<index_t> body_edge_starts;
    std::vector<index_t> body_face_starts;
};

} // namespace ipc::rigid