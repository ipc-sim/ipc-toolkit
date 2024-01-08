#pragma once

#include <ipc/collision_mesh.hpp>
#include "smooth_collisions.hpp"

#include <tbb/enumerable_thread_specific.h>

#include <Eigen/Core>

namespace ipc {

class SmoothCollisionsBuilder {
public:
    SmoothCollisionsBuilder() = default;

    void add_edge_vertex_collisions(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const std::vector<EdgeVertexCandidate>& candidates,
        const std::function<bool(double)>& is_active,
        const size_t start_i,
        const size_t end_i);

    // ------------------------------------------------------------------------

    static void merge(
        const tbb::enumerable_thread_specific<SmoothCollisionsBuilder>& local_storage,
        SmoothCollisions& merged_collisions);

    // -------------------------------------------------------------------------

    static void add_edge_vertex_collision(
        const SmoothEdgeVertexCollision& ev_collision,
        unordered_map<SmoothEdgeVertexCollision, long>& ev_to_id,
        std::vector<SmoothEdgeVertexCollision>& ev_collisions);

    void add_edge_vertex_collision(
        const long edge_id,
        const long vertex_id,
        const double weight,
        const Eigen::SparseVector<double>& weight_gradient,
        const double local_eps)
    {
        add_edge_vertex_collision(
            SmoothEdgeVertexCollision(edge_id, vertex_id, weight, weight_gradient, local_eps),
            ev_to_id, ev_collisions);
    }

    void add_edge_vertex_collision(
        const CollisionMesh& mesh,
        const EdgeVertexCandidate& candidate,
        const PointEdgeDistanceType dtype,
        const double weight,
        const Eigen::SparseVector<double>& weight_gradient);

    void add_edge_edge_collisions(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const std::vector<EdgeEdgeCandidate>& candidates,
        const std::function<bool(double)>& is_active,
        const size_t start_i,
        const size_t end_i,
        const SurfaceQuadratureType quad_type);

    static void add_edge_edge_collision(
        const SmoothEdgeEdgeCollision& ee_collision,
        unordered_map<SmoothEdgeEdgeCollision, long>& ee_to_id_,
        std::vector<SmoothEdgeEdgeCollision>& ee_collisions_);

    // -------------------------------------------------------------------------

    // Store the indices to pairs to avoid duplicates.
    // unordered_map<VertexVertexCollision, long> vv_to_id;
    unordered_map<SmoothEdgeVertexCollision, long> ev_to_id;
    unordered_map<SmoothEdgeEdgeCollision, long> ee_to_id;

    // Constructed collisions
    // std::vector<VertexVertexCollision> vv_collisions;
    std::vector<SmoothEdgeVertexCollision> ev_collisions;
    std::vector<SmoothEdgeEdgeCollision> ee_collisions;
    // std::vector<FaceVertexCollision> fv_collisions;
    // std::vector<PlaneVertexCollision> pv_collisions;
};

} // namespace ipc