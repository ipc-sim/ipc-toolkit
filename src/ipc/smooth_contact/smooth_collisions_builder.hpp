#pragma once

#include <ipc/collision_mesh.hpp>
#include "smooth_collisions.hpp"

#include <tbb/enumerable_thread_specific.h>

#include <Eigen/Core>

namespace ipc {

template <int dim, class TCollision>
class SmoothCollisionsBuilder {
public:
    SmoothCollisionsBuilder() = default;

    // only for 2D, transform edge-vertex to edge-edge
    void add_edge_vertex_collisions(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const std::vector<EdgeVertexCandidate>& candidates,
        const std::function<bool(double)>& is_active,
        const size_t start_i,
        const size_t end_i,
        const double dhat,
        const bool use_adaptive_eps);

    void add_neighbor_edge_collisions(
        const CollisionMesh& mesh,
        const size_t start_i,
        const size_t end_i,
        const bool use_adaptive_eps);

    // ------------------------------------------------------------------------

    static void merge(
        const tbb::enumerable_thread_specific<SmoothCollisionsBuilder<dim, TCollision>>& local_storage,
        SmoothCollisions<dim>& merged_collisions);

    // -------------------------------------------------------------------------

    // only for 3D, transform edge-edge to face-face
    void add_edge_edge_collisions(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const std::vector<EdgeEdgeCandidate>& candidates,
        const std::function<bool(double)>& is_active,
        const size_t start_i,
        const size_t end_i,
        const double dhat,
        const bool use_adaptive_eps = false);

    // only for 3D, transform face-vertex to face-face
    void add_face_vertex_collisions(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const std::vector<FaceVertexCandidate>& candidates,
        const std::function<bool(double)>& is_active,
        const size_t start_i,
        const size_t end_i);

    static void add_collision(
        const TCollision& collision,
        unordered_map<TCollision, long>& cc_to_id_,
        std::vector<std::shared_ptr<TCollision>>& collisions_);

    // -------------------------------------------------------------------------

    // Store the indices to pairs to avoid duplicates.
    unordered_map<TCollision, long> cc_to_id;

    // Constructed collisions
    std::vector<std::shared_ptr<TCollision>> collisions;
};

} // namespace ipc