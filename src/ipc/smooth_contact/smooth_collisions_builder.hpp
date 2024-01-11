#pragma once

#include <ipc/collision_mesh.hpp>
#include "smooth_collisions.hpp"
#include <ipc/utils/unordered_tuple.hpp>
#include <tbb/enumerable_thread_specific.h>

#include <Eigen/Core>

namespace ipc {

template <int dim>
class SmoothCollisionsBuilder {
    using TCollision = typename SmoothCollision<dim>::Type;
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
        const bool use_adaptive_dhat);

    void add_neighbor_edge_collisions(
        const CollisionMesh& mesh,
        const size_t start_i,
        const size_t end_i,
        const bool use_adaptive_dhat);

    // ------------------------------------------------------------------------

    static void merge(
        const CollisionMesh &mesh,
        const tbb::enumerable_thread_specific<SmoothCollisionsBuilder<dim>>& local_storage,
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
        const bool use_adaptive_dhat = false);

    // only for 3D, transform face-vertex to face-face
    void add_face_vertex_collisions(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const std::vector<FaceVertexCandidate>& candidates,
        const std::function<bool(double)>& is_active,
        const size_t start_i,
        const size_t end_i);

    static void add_collision(
        const CollisionMesh &mesh,
        const unordered_tuple& pair,
        unordered_map<unordered_tuple, long>& cc_to_id_,
        std::vector<TCollision>& collisions_);

    // -------------------------------------------------------------------------

    // Store the indices to pairs to avoid duplicates.
    unordered_map<unordered_tuple, long> cc_to_id;

    // Constructed collisions
    std::vector<TCollision> collisions;
};

} // namespace ipc