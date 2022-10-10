#pragma once

#include <ipc/collision_mesh.hpp>
#include <ipc/broad_phase/collision_candidate.hpp>
#include <ipc/collisions/constraints.hpp>

#include <tbb/enumerable_thread_specific.h>

#include <Eigen/Core>

namespace ipc {

class CollisionConstraintsBuilder {
public:
    CollisionConstraintsBuilder(
        const bool use_convergent_formulation,
        const bool compute_shape_derivatives);

    void add_edge_vertex_constraints(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& V,
        const std::vector<EdgeVertexCandidate>& candidates,
        const std::function<bool(double)>& is_active,
        const size_t start_i,
        const size_t end_i);

    void add_edge_edge_constraints(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& V,
        const std::vector<EdgeEdgeCandidate>& candidates,
        const std::function<bool(double)>& is_active,
        const size_t start_i,
        const size_t end_i);

    void add_face_vertex_constraints(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& V,
        const std::vector<FaceVertexCandidate>& candidates,
        const std::function<bool(double)>& is_active,
        const size_t start_i,
        const size_t end_i);

    static void merge(
        const tbb::enumerable_thread_specific<CollisionConstraintsBuilder>&
            local_storage,
        CollisionConstraints& constraints);

protected:
    static void add_vertex_vertex_constraint(
        const long v0i,
        const long v1i,
        const double weight,
        const Eigen::SparseVector<double>& weight_gradient,
        unordered_map<VertexVertexConstraint, long>& vv_to_index,
        std::vector<VertexVertexConstraint>& vv_constraints);

    static void add_edge_vertex_constraint(
        const long ei,
        const long vi,
        const double weight,
        const Eigen::SparseVector<double>& weight_gradient,
        unordered_map<EdgeVertexConstraint, long>& ev_to_index,
        std::vector<EdgeVertexConstraint>& ev_constraints);

    void add_vertex_vertex_constraint(
        const long v0i,
        const long v1i,
        const double weight,
        const Eigen::SparseVector<double>& weight_gradient)
    {
        add_vertex_vertex_constraint(
            v0i, v1i, weight, weight_gradient, vv_to_index,
            constraints.vv_constraints);
    }

    void add_edge_vertex_constraint(
        const long ei,
        const long vi,
        const double weight,
        const Eigen::SparseVector<double>& weight_gradient)
    {
        add_edge_vertex_constraint(
            ei, vi, weight, weight_gradient, ev_to_index,
            constraints.ev_constraints);
    }

    // Store the indices to VV and EV pairs to avoid duplicates.
    unordered_map<VertexVertexConstraint, long> vv_to_index;
    unordered_map<EdgeVertexConstraint, long> ev_to_index;
    CollisionConstraints constraints;
    bool use_convergent_formulation;
    bool compute_shape_derivatives;
};

} // namespace ipc