// NOTE: This is an internal header file, not meant to be used outside of the
// IPC Toolkit library. It includes TBB which is a private dependency of the IPC
// Toolkit library. To use this outside of the library, one needs to link
// against TBB::tbb.

#pragma once

#include <ipc/collision_mesh.hpp>
#include <ipc/collisions/normal/normal_collisions.hpp>
#include <ipc/utils/unordered_map_and_set.hpp>

#include <Eigen/Core>
#include <tbb/enumerable_thread_specific.h>

namespace ipc {

class NormalCollisionsBuilder {
public:
    NormalCollisionsBuilder(
        const bool use_area_weighting,
        const bool enable_shape_derivatives,
        const bool use_ogc);

    void add_vertex_vertex_collision(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        const VertexVertexCandidate& candidate,
        const std::function<bool(double)>& is_active);

    void add_edge_vertex_collision(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        const EdgeVertexCandidate& candidate,
        const std::function<bool(double)>& is_active);

    void add_edge_edge_collision(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        const EdgeEdgeCandidate& candidate,
        const std::function<bool(double)>& is_active);

    void add_face_vertex_collision(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        const FaceVertexCandidate& candidate,
        const std::function<bool(double)>& is_active);

    // ------------------------------------------------------------------------
    // Duplicate removal functions

    void add_edge_vertex_negative_vertex_vertex_collision(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        const VertexVertexCandidate& candidate);

    void add_face_vertex_positive_vertex_vertex_collision(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        const VertexVertexCandidate& candidate);

    void add_face_vertex_negative_edge_vertex_collision(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        const EdgeVertexCandidate& candidate);

    void add_edge_edge_negative_edge_vertex_collision(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        const EdgeVertexCandidate& candidate);

    // ------------------------------------------------------------------------

    static void merge(
        const tbb::enumerable_thread_specific<NormalCollisionsBuilder>&
            local_storage,
        NormalCollisions& merged_collisions);

    // -------------------------------------------------------------------------
protected:
    static void add_vertex_vertex_collision(
        const VertexVertexNormalCollision& vv_collision,
        unordered_map<VertexVertexNormalCollision, size_t>& vv_to_id,
        std::vector<VertexVertexNormalCollision>& vv_collisions);

    void add_vertex_vertex_collision(
        const index_t vertex0_id,
        const index_t vertex1_id,
        const double weight,
        const Eigen::SparseVector<double>& weight_gradient)
    {
        add_vertex_vertex_collision(
            VertexVertexNormalCollision(
                vertex0_id, vertex1_id, weight, weight_gradient),
            vv_to_id, vv_collisions);
    }

    // -------------------------------------------------------------------------

    static void add_edge_vertex_collision(
        const EdgeVertexNormalCollision& ev_collision,
        unordered_map<EdgeVertexNormalCollision, size_t>& ev_to_id,
        std::vector<EdgeVertexNormalCollision>& ev_collisions);

    void add_edge_vertex_collision(
        const index_t edge_id,
        const index_t vertex_id,
        const double weight,
        const Eigen::SparseVector<double>& weight_gradient)
    {
        add_edge_vertex_collision(
            EdgeVertexNormalCollision(
                edge_id, vertex_id, weight, weight_gradient),
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
        const EdgeEdgeNormalCollision& ee_collision,
        unordered_map<EdgeEdgeNormalCollision, size_t>& ee_to_id,
        std::vector<EdgeEdgeNormalCollision>& ee_collisions);

    void add_edge_edge_collision(
        const index_t edge0_id,
        const index_t edge1_id,
        const double eps_x,
        const double weight,
        const Eigen::SparseVector<double>& weight_gradient,
        const EdgeEdgeDistanceType dtype)
    {
        add_edge_edge_collision(
            EdgeEdgeNormalCollision(
                edge0_id, edge1_id, eps_x, weight, weight_gradient, dtype),
            ee_to_id, ee_collisions);
    }

    // -------------------------------------------------------------------------

    // Store the indices to pairs to avoid duplicates.
    unordered_map<VertexVertexNormalCollision, size_t> vv_to_id;
    unordered_map<EdgeVertexNormalCollision, size_t> ev_to_id;
    unordered_map<EdgeEdgeNormalCollision, size_t> ee_to_id;

    // Constructed collisions
    std::vector<VertexVertexNormalCollision> vv_collisions;
    std::vector<EdgeVertexNormalCollision> ev_collisions;
    std::vector<EdgeEdgeNormalCollision> ee_collisions;
    std::vector<FaceVertexNormalCollision> fv_collisions;
    // std::vector<PlaneVertexNormalCollision> pv_collisions;

    const bool use_area_weighting;
    const bool enable_shape_derivatives;
    const bool use_ogc;
};

} // namespace ipc