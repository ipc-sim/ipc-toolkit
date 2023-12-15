#pragma once

#include <ipc/collision_mesh.hpp>
#include <ipc/collisions/collisions.hpp>

#include <tbb/enumerable_thread_specific.h>

#include <Eigen/Core>

namespace ipc {

class CollisionsBuilder {
public:
    CollisionsBuilder(
        const bool use_convergent_formulation,
        const bool are_shape_derivatives_enabled);

    void add_vertex_vertex_collisions(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const std::vector<VertexVertexCandidate>& candidates,
        const std::function<bool(double)>& is_active,
        const size_t start_i,
        const size_t end_i);

    void add_edge_vertex_collisions(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const std::vector<EdgeVertexCandidate>& candidates,
        const std::function<bool(double)>& is_active,
        const size_t start_i,
        const size_t end_i);

    void add_edge_edge_collisions(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const std::vector<EdgeEdgeCandidate>& candidates,
        const std::function<bool(double)>& is_active,
        const size_t start_i,
        const size_t end_i);

    void add_face_vertex_collisions(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const std::vector<FaceVertexCandidate>& candidates,
        const std::function<bool(double)>& is_active,
        const size_t start_i,
        const size_t end_i);

    // ------------------------------------------------------------------------
    // Duplicate removal functions

    void add_edge_vertex_negative_vertex_vertex_collisions(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const std::vector<VertexVertexCandidate>& candidates,
        const size_t start_i,
        const size_t end_i);

    void add_face_vertex_positive_vertex_vertex_collisions(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const std::vector<VertexVertexCandidate>& candidates,
        const size_t start_i,
        const size_t end_i);

    void add_face_vertex_negative_edge_vertex_collisions(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const std::vector<EdgeVertexCandidate>& candidates,
        const size_t start_i,
        const size_t end_i);

    void add_edge_edge_negative_edge_vertex_collisions(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const std::vector<EdgeVertexCandidate>& candidates,
        const size_t start_i,
        const size_t end_i);

    // ------------------------------------------------------------------------

    static void merge(
        const tbb::enumerable_thread_specific<CollisionsBuilder>& local_storage,
        Collisions& merged_collisions);

    // -------------------------------------------------------------------------
protected:
    static void add_vertex_vertex_collision(
        const VertexVertexCollision& vv_collision,
        unordered_map<VertexVertexCollision, long>& vv_to_id,
        std::vector<VertexVertexCollision>& vv_collisions);

    void add_vertex_vertex_collision(
        const long vertex0_id,
        const long vertex1_id,
        const double weight,
        const Eigen::SparseVector<double>& weight_gradient)
    {
        add_vertex_vertex_collision(
            VertexVertexCollision(
                vertex0_id, vertex1_id, weight, weight_gradient),
            vv_to_id, vv_collisions);
    }

    // -------------------------------------------------------------------------

    static void add_edge_vertex_collision(
        const EdgeVertexCollision& ev_collision,
        unordered_map<EdgeVertexCollision, long>& ev_to_id,
        std::vector<EdgeVertexCollision>& ev_collisions);

    void add_edge_vertex_collision(
        const long edge_id,
        const long vertex_id,
        const double weight,
        const Eigen::SparseVector<double>& weight_gradient)
    {
        add_edge_vertex_collision(
            EdgeVertexCollision(edge_id, vertex_id, weight, weight_gradient),
            ev_to_id, ev_collisions);
    }

    void add_edge_vertex_collision(
        const CollisionMesh& mesh,
        const EdgeVertexCandidate& candidate,
        const PointEdgeDistanceType dtype,
        const double weight,
        const Eigen::SparseVector<double>& weight_gradient);

    // -------------------------------------------------------------------------

    static void add_edge_edge_collision(
        const EdgeEdgeCollision& ee_collision,
        unordered_map<EdgeEdgeCollision, long>& ee_to_id,
        std::vector<EdgeEdgeCollision>& ee_collisions);

    void add_edge_edge_collision(
        const long edge0_id,
        const long edge1_id,
        const double eps_x,
        const double weight,
        const Eigen::SparseVector<double>& weight_gradient,
        const EdgeEdgeDistanceType dtype)
    {
        add_edge_edge_collision(
            EdgeEdgeCollision(
                edge0_id, edge1_id, eps_x, weight, weight_gradient, dtype),
            ee_to_id, ee_collisions);
    }

    // -------------------------------------------------------------------------

    // Store the indices to pairs to avoid duplicates.
    unordered_map<VertexVertexCollision, long> vv_to_id;
    unordered_map<EdgeVertexCollision, long> ev_to_id;
    unordered_map<EdgeEdgeCollision, long> ee_to_id;

    // Constructed collisions
    std::vector<VertexVertexCollision> vv_collisions;
    std::vector<EdgeVertexCollision> ev_collisions;
    std::vector<EdgeEdgeCollision> ee_collisions;
    std::vector<FaceVertexCollision> fv_collisions;
    // std::vector<PlaneVertexCollision> pv_collisions;

    const bool use_convergent_formulation;
    const bool should_compute_weight_gradient;
};

} // namespace ipc