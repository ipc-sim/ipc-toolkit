#pragma once

#include <ipc/collision_mesh.hpp>
#include <ipc/collisions/collision_constraints.hpp>

#include <tbb/enumerable_thread_specific.h>

#include <Eigen/Core>

namespace ipc {

class CollisionConstraintsBuilder {
public:
    CollisionConstraintsBuilder(
        const bool use_convergent_formulation,
        const bool are_shape_derivatives_enabled);

    void add_vertex_vertex_constraints(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const std::vector<VertexVertexCandidate>& candidates,
        const std::function<bool(double)>& is_active,
        const size_t start_i,
        const size_t end_i);

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

    // ------------------------------------------------------------------------
    // Duplicate removal functions

    void add_edge_vertex_negative_vertex_vertex_constraints(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const std::vector<VertexVertexCandidate>& candidates,
        const size_t start_i,
        const size_t end_i);

    void add_face_vertex_positive_vertex_vertex_constraints(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const std::vector<VertexVertexCandidate>& candidates,
        const size_t start_i,
        const size_t end_i);

    void add_face_vertex_negative_edge_vertex_constraints(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const std::vector<EdgeVertexCandidate>& candidates,
        const size_t start_i,
        const size_t end_i);

    void add_edge_edge_negative_edge_vertex_constraints(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const std::vector<EdgeVertexCandidate>& candidates,
        const size_t start_i,
        const size_t end_i);

    // ------------------------------------------------------------------------

    static void merge(
        const tbb::enumerable_thread_specific<CollisionConstraintsBuilder>&
            local_storage,
        CollisionConstraints& merged_constraints);

    // -------------------------------------------------------------------------
protected:
    static void add_vertex_vertex_constraint(
        const VertexVertexConstraint& vv_constraint,
        unordered_map<VertexVertexConstraint, long>& vv_to_id,
        std::vector<VertexVertexConstraint>& vv_constraints);

    void add_vertex_vertex_constraint(
        const long vertex0_id,
        const long vertex1_id,
        const double weight,
        const Eigen::SparseVector<double>& weight_gradient)
    {
        add_vertex_vertex_constraint(
            VertexVertexConstraint(
                vertex0_id, vertex1_id, weight, weight_gradient),
            vv_to_id, vv_constraints);
    }

    // -------------------------------------------------------------------------

    static void add_edge_vertex_constraint(
        const EdgeVertexConstraint& ev_constraint,
        unordered_map<EdgeVertexConstraint, long>& ev_to_id,
        std::vector<EdgeVertexConstraint>& ev_constraints);

    void add_edge_vertex_constraint(
        const long edge_id,
        const long vertex_id,
        const double weight,
        const Eigen::SparseVector<double>& weight_gradient)
    {
        add_edge_vertex_constraint(
            EdgeVertexConstraint(edge_id, vertex_id, weight, weight_gradient),
            ev_to_id, ev_constraints);
    }

    void add_edge_vertex_constraint(
        const CollisionMesh& mesh,
        const EdgeVertexCandidate& candidate,
        const PointEdgeDistanceType dtype,
        const double weight,
        const Eigen::SparseVector<double>& weight_gradient);

    // -------------------------------------------------------------------------

    static void add_edge_edge_constraint(
        const EdgeEdgeConstraint& ee_constraint,
        unordered_map<EdgeEdgeConstraint, long>& ee_to_id,
        std::vector<EdgeEdgeConstraint>& ee_constraints);

    void add_edge_edge_constraint(
        const long edge0_id,
        const long edge1_id,
        const double eps_x,
        const double weight,
        const Eigen::SparseVector<double>& weight_gradient,
        const EdgeEdgeDistanceType dtype)
    {
        add_edge_edge_constraint(
            EdgeEdgeConstraint(
                edge0_id, edge1_id, eps_x, weight, weight_gradient, dtype),
            ee_to_id, ee_constraints);
    }

    // -------------------------------------------------------------------------

    // Store the indices to pairs to avoid duplicates.
    unordered_map<VertexVertexConstraint, long> vv_to_id;
    unordered_map<EdgeVertexConstraint, long> ev_to_id;
    unordered_map<EdgeEdgeConstraint, long> ee_to_id;

    // Constructed constraints
    std::vector<VertexVertexConstraint> vv_constraints;
    std::vector<EdgeVertexConstraint> ev_constraints;
    std::vector<EdgeEdgeConstraint> ee_constraints;
    std::vector<FaceVertexConstraint> fv_constraints;
    // std::vector<PlaneVertexConstraint> pv_constraints;

    const bool use_convergent_formulation;
    const bool should_compute_weight_gradient;
};

} // namespace ipc