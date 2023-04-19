#pragma once

#include <ipc/collision_mesh.hpp>
#include <ipc/candidates/vertex_vertex.hpp>
#include <ipc/candidates/edge_vertex.hpp>
#include <ipc/collisions/collision_constraints.hpp>

#include <tbb/enumerable_thread_specific.h>

#include <Eigen/Core>

namespace ipc {

class CollisionConstraintsBuilder {
public:
    CollisionConstraintsBuilder(const CollisionConstraints& empty_constraints);

    void add_edge_vertex_constraints(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const std::vector<EdgeVertexCandidate>& candidates,
        const std::function<bool(double)>& is_active,
        const size_t start_i,
        const size_t end_i);

    void add_edge_edge_constraints(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const std::vector<EdgeEdgeCandidate>& candidates,
        const std::function<bool(double)>& is_active,
        const size_t start_i,
        const size_t end_i);

    void add_face_vertex_constraints(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const std::vector<FaceVertexCandidate>& candidates,
        const std::function<bool(double)>& is_active,
        const size_t start_i,
        const size_t end_i);

    static void merge(
        const tbb::enumerable_thread_specific<CollisionConstraintsBuilder>&
            local_storage,
        CollisionConstraints& merged_constraints);

protected:
    static void add_vertex_vertex_constraint(
        const long v0i,
        const long v1i,
        const double weight,
        const Eigen::SparseVector<double>& weight_gradient,
        unordered_map<VertexVertexConstraint, long>& vv_to_id,
        std::vector<VertexVertexConstraint>& vv_constraints);

    static void add_edge_vertex_constraint(
        const long ei,
        const long vi,
        const double weight,
        const Eigen::SparseVector<double>& weight_gradient,
        unordered_map<EdgeVertexConstraint, long>& ev_to_id,
        std::vector<EdgeVertexConstraint>& ev_constraints);

    void add_vertex_vertex_constraint(
        const long v0i,
        const long v1i,
        const double weight,
        const Eigen::SparseVector<double>& weight_gradient)
    {
        add_vertex_vertex_constraint(
            v0i, v1i, weight, weight_gradient, vv_to_id,
            constraints.vv_constraints);
    }

    void add_edge_vertex_constraint(
        const long ei,
        const long vi,
        const double weight,
        const Eigen::SparseVector<double>& weight_gradient)
    {
        add_edge_vertex_constraint(
            ei, vi, weight, weight_gradient, ev_to_id,
            constraints.ev_constraints);
    }

    bool use_convergent_formulation() const
    {
        return constraints.use_convergent_formulation();
    }

    bool should_compute_weight_gradient() const
    {
        return constraints.are_shape_derivatives_enabled();
    }

    // Store the indices to VV and EV pairs to avoid duplicates.
    unordered_map<VertexVertexConstraint, long> vv_to_id;
    unordered_map<EdgeVertexConstraint, long> ev_to_id;
    CollisionConstraints constraints;
};

} // namespace ipc