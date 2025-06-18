#include "rigid_bodies.hpp"

#include <ipc/utils/logger.hpp>

namespace ipc::rigid {

RigidBodies RigidBodies::build_from_meshes(
    const std::vector<Eigen::MatrixXd>& rest_positions,
    const std::vector<Eigen::MatrixXi>& edges,
    const std::vector<Eigen::MatrixXi>& faces,
    const std::vector<double>& densities,
    std::vector<Pose>& initial_poses)
{
    assert(rest_positions.size() == edges.size());
    assert(rest_positions.size() == faces.size());
    assert(rest_positions.size() == densities.size());

    size_t num_vertices = 0, num_edges = 0, num_faces = 0;
    std::vector<index_t> body_vertex_starts(rest_positions.size() + 1);
    std::vector<index_t> body_edge_starts(edges.size() + 1);
    std::vector<index_t> body_face_starts(faces.size() + 1);
    body_vertex_starts[0] = body_edge_starts[0] = body_face_starts[0] = 0;

    for (size_t i = 0; i < rest_positions.size(); ++i) {
        body_vertex_starts[i + 1] = (num_vertices += rest_positions[i].rows());
        body_edge_starts[i + 1] = (num_edges += edges[i].rows());
        body_face_starts[i + 1] = (num_faces += faces[i].rows());
    }

    Eigen::MatrixXd concat_rest_positions(
        num_vertices, rest_positions[0].cols());
    Eigen::MatrixXi concat_edges(num_edges, 2);
    Eigen::MatrixXi concat_faces(num_faces, 3);

    for (size_t i = 0; i < rest_positions.size(); ++i) {
        concat_rest_positions.middleRows(
            body_vertex_starts[i], rest_positions[i].rows()) =
            rest_positions[i];
        concat_edges.middleRows(body_edge_starts[i], edges[i].rows()) =
            edges[i].array() + body_vertex_starts[i];
        concat_faces.middleRows(body_face_starts[i], faces[i].rows()) =
            faces[i].array() + body_vertex_starts[i];
    }

    return RigidBodies(
        concat_rest_positions, concat_edges, concat_faces, body_vertex_starts,
        body_edge_starts, body_face_starts, densities, initial_poses);
}

RigidBodies::RigidBodies(
    Eigen::ConstRef<Eigen::MatrixXd> _rest_positions,
    Eigen::ConstRef<Eigen::MatrixXi> _edges,
    Eigen::ConstRef<Eigen::MatrixXi> _faces,
    const std::vector<index_t>& _body_vertex_starts,
    const std::vector<index_t>& _body_edge_starts,
    const std::vector<index_t>& _body_face_starts,
    const std::vector<double>& densities,
    std::vector<Pose>& initial_poses)
    : CollisionMesh(_rest_positions, _edges, _faces)
    , body_vertex_starts(std::move(_body_vertex_starts))
    , body_edge_starts(std::move(_body_edge_starts))
    , body_face_starts(std::move(_body_face_starts))
{
    assert(body_vertex_starts.size() == body_edge_starts.size());
    assert(body_edge_starts.size() == body_face_starts.size());
    assert(body_vertex_starts.back() == num_vertices());
    assert(body_edge_starts.back() == num_edges());
    assert(body_face_starts.back() == num_faces());
    assert(initial_poses.size() == body_vertex_starts.size() - 1);

    bodies.reserve(body_vertex_starts.size() - 1);
    for (size_t i = 0; i < body_vertex_starts.size() - 1; ++i) {
        bodies.emplace_back(
            m_rest_positions.middleRows(
                body_vertex_starts[i],
                body_vertex_starts[i + 1] - body_vertex_starts[i]),
            edges().middleRows(
                body_edge_starts[i],
                body_edge_starts[i + 1] - body_edge_starts[i]),
            faces().middleRows(
                body_face_starts[i],
                body_face_starts[i + 1] - body_face_starts[i]),
            densities[i], initial_poses[i]);
        logger().info(
            "Initial pose: position={}, rotation={}",
            initial_poses[i].position.transpose(),
            initial_poses[i].rotation.transpose());
    }
}

} // namespace ipc::rigid