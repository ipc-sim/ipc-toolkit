#pragma once

#include <Eigen/Core>

#include <ipc/broad_phase/collision_candidate.hpp>

namespace ipc {

/// Find all static collisions in one time step using brute-force.
void detect_collision_candidates_brute_force(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    Candidates& candidates,
    bool detect_edge_vertex = false,
    bool detect_edge_edge = true,
    bool detect_face_vertex = true,
    bool perform_aabb_check = false,
    double aabb_inflation_radius = 0,
    bool ignore_codimensional_vertices = false,
    const std::function<bool(size_t, size_t)>& can_collide =
        [](size_t, size_t) { return true; });

void detect_edge_vertex_collision_candidates_brute_force(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    std::vector<EdgeVertexCandidate>& ev_candidates,
    bool perform_aabb_check = false,
    double aabb_inflation_radius = 0,
    const std::function<bool(size_t, size_t)>& can_collide =
        [](size_t, size_t) { return true; });

void detect_edge_edge_collision_candidates_brute_force(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    std::vector<EdgeEdgeCandidate>& ee_candidates,
    bool perform_aabb_check = false,
    double aabb_inflation_radius = 0,
    const std::function<bool(size_t, size_t)>& can_collide =
        [](size_t, size_t) { return true; });

void detect_face_vertex_collision_candidates_brute_force(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    std::vector<FaceVertexCandidate>& fv_candidates,
    bool perform_aabb_check = false,
    double aabb_inflation_radius = 0,
    const std::function<bool(size_t, size_t)>& can_collide =
        [](size_t, size_t) { return true; });

void detect_edge_face_collision_candidates_brute_force(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    std::vector<EdgeFaceCandidate>& ef_candidates,
    bool perform_aabb_check = false,
    double aabb_inflation_radius = 0,
    const std::function<bool(size_t, size_t)>& can_collide =
        [](size_t, size_t) { return true; });

/// Find all continous collisions in one time step using brute-force.
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
    bool ignore_codimensional_vertices = false,
    const std::function<bool(size_t, size_t)>& can_collide =
        [](size_t, size_t) { return true; });

void detect_edge_vertex_collision_candidates_brute_force(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    std::vector<EdgeVertexCandidate>& ev_candidates,
    bool perform_aabb_check = false,
    double aabb_inflation_radius = 0,
    const std::function<bool(size_t, size_t)>& can_collide =
        [](size_t, size_t) { return true; });

void detect_edge_edge_collision_candidates_brute_force(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    std::vector<EdgeEdgeCandidate>& ee_candidates,
    bool perform_aabb_check = false,
    double aabb_inflation_radius = 0,
    const std::function<bool(size_t, size_t)>& can_collide =
        [](size_t, size_t) { return true; });

void detect_face_vertex_collision_candidates_brute_force(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& F,
    std::vector<FaceVertexCandidate>& fv_candidates,
    bool perform_aabb_check = false,
    double aabb_inflation_radius = 0,
    const std::function<bool(size_t, size_t)>& can_collide =
        [](size_t, size_t) { return true; });

} // namespace ipc
