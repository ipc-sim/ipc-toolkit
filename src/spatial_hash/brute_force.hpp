#pragma once

#include <Eigen/Core>

#include <ipc/spatial_hash/collision_candidate.hpp>

namespace ipc {

// Find all edge-vertex collisions in one time step using brute-force
// comparisons of all E and all vertices.
void detect_collision_candidates_brute_force(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    Candidates& candidates,
    bool detect_edge_vertex = false,
    bool detect_edge_edge = true,
    bool detect_face_vertex = true,
    bool perform_aabb_check = false,
    double aabb_inflation_radius = 0,
    const Eigen::VectorXi& group_ids = Eigen::VectorXi());

void detect_edge_vertex_collision_candidates_brute_force(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    std::vector<EdgeVertexCandidate>& ev_candidates,
    bool perform_aabb_check = false,
    double aabb_inflation_radius = 0,
    const Eigen::VectorXi& group_ids = Eigen::VectorXi());

void detect_edge_edge_collision_candidates_brute_force(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    std::vector<EdgeEdgeCandidate>& ee_candidates,
    bool perform_aabb_check = false,
    double aabb_inflation_radius = 0,
    const Eigen::VectorXi& group_ids = Eigen::VectorXi());

void detect_face_vertex_collision_candidates_brute_force(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& F,
    std::vector<FaceVertexCandidate>& fv_candidates,
    bool perform_aabb_check = false,
    double aabb_inflation_radius = 0,
    const Eigen::VectorXi& group_ids = Eigen::VectorXi());

} // namespace ipc
