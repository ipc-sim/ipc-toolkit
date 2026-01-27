#pragma once

#include <ipc/collision_mesh.hpp>
#include <ipc/candidates/edge_edge.hpp>
#include <ipc/candidates/edge_vertex.hpp>
#include <ipc/candidates/face_vertex.hpp>

namespace ipc::ogc {

/// @brief Check if point `x` is in the feasible region of vertex `vi`.
/// @param mesh Collision mesh containing the vertex adjacencies
/// @param vertices Matrix of current vertex positions (rowwise)
/// @param x Position of the point to check
/// @param vi Index of the vertex for which to check the feasible region
/// @return True if the point is in the feasible region, false otherwise
bool check_vertex_feasible_region(
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices,
    Eigen::ConstRef<VectorMax3d> x,
    const index_t vi);

/// @brief Check if vertex `xi` is in the feasible region of vertex `vi`.
/// @param mesh Collision mesh containing the vertex adjacencies
/// @param vertices Matrix of current vertex positions (rowwise)
/// @param xi Index of the vertex to check
/// @param vi Index of the vertex for which to check the feasible region
/// @return True if the vertex is in the feasible region, false otherwise
inline bool check_vertex_feasible_region(
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices,
    const index_t xi,
    const index_t vi)
{
    return check_vertex_feasible_region(mesh, vertices, vertices.row(xi), vi);
}

/// @brief Check if vertex `xi` is in the feasible region of edge `ei`.
/// @param mesh Collision mesh containing the edge adjacencies
/// @param vertices Matrix of current vertex positions (rowwise)
/// @param xi Index of the vertex to check
/// @param ei Index of the edge for which to check the feasible region
/// @return True if the vertex is in the feasible region, false otherwise
bool check_edge_feasible_region(
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices,
    const index_t xi,
    const index_t ei);

/// @brief Check if the edge-vertex candidate is feasible.
/// @param mesh Collision mesh containing the edge adjacencies
/// @param vertices Matrix of current vertex positions (rowwise)
/// @param candidate Edge-vertex candidate to check
/// @param dtype Edge-vertex distance type to use for the check.
/// @return True if the edge-vertex candidate is in the feasible region, false otherwise
bool is_edge_vertex_feasible(
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices,
    const EdgeVertexCandidate& candidate,
    PointEdgeDistanceType dtype = PointEdgeDistanceType::AUTO);

/// @brief Check if the edge-edge candidate is feasible.
/// @param mesh Collision mesh containing the edge adjacencies
/// @param vertices Matrix of current vertex positions (rowwise)
/// @param candidate Edge-edge candidate to check
/// @param dtype Edge-edge distance type to use for the check.
/// @return True if the edge-edge candidate is in the feasible region, false otherwise
bool is_edge_edge_feasible(
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices,
    const EdgeEdgeCandidate& candidate,
    EdgeEdgeDistanceType dtype = EdgeEdgeDistanceType::AUTO);

/// @brief Check if the face-vertex candidate is feasible.
/// @param mesh Collision mesh containing the face adjacencies
/// @param vertices Matrix of current vertex positions (rowwise)
/// @param candidate Face-vertex candidate to check
/// @param dtype Point-triangle distance type to use for the check.
/// @return True if the face-vertex candidate is in the feasible region, false otherwise
bool is_face_vertex_feasible(
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices,
    const FaceVertexCandidate& candidate,
    PointTriangleDistanceType dtype = PointTriangleDistanceType::AUTO);

} // namespace ipc::ogc