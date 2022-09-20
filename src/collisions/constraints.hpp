#pragma once

#include <ipc/collision_mesh.hpp>
#include <ipc/collisions/collision_constraint.hpp>
#include <ipc/collisions/vertex_vertex.hpp>
#include <ipc/collisions/edge_vertex.hpp>
#include <ipc/collisions/edge_edge.hpp>
#include <ipc/collisions/face_vertex.hpp>
#include <ipc/collisions/plane_vertex.hpp>
#include <ipc/broad_phase/broad_phase.hpp>

#include <Eigen/Core>

#include <tbb/enumerable_thread_specific.h>

#include <vector>

namespace ipc {

struct Constraints {
    std::vector<VertexVertexConstraint> vv_constraints;
    std::vector<EdgeVertexConstraint> ev_constraints;
    std::vector<EdgeEdgeConstraint> ee_constraints;
    std::vector<FaceVertexConstraint> fv_constraints;
    std::vector<PlaneVertexConstraint> pv_constraints;
    bool use_convergent_formulation = false;
    bool compute_shape_derivatives = false;

    /// @brief Construct a set of constraints used to compute the barrier potential.
    /// @param mesh The collision mesh.
    /// @param V Vertices of the collision mesh.
    /// @param dhat The activation distance of the barrier.
    /// @param dmin Minimum distance.
    /// @param method Broad-phase method to use.
    void build(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& V,
        const double dhat,
        const double dmin = 0,
        const BroadPhaseMethod method = BroadPhaseMethod::HASH_GRID);

    /// @brief Construct a set of constraints used to compute the barrier potential.
    /// @param candidates Distance candidates from which the constraint set is built.
    /// @param mesh The collision mesh.
    /// @param V Vertices of the collision mesh.
    /// @param dhat The activation distance of the barrier.
    /// @param  dmin  Minimum distance.
    void build(
        const Candidates& candidates,
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& V,
        const double dhat,
        const double dmin = 0);

    size_t size() const;

    bool empty() const;

    void clear();

    CollisionConstraint& operator[](size_t idx);
    const CollisionConstraint& operator[](size_t idx) const;

protected:
    struct Builder;

    void edge_vertex_candiates_to_constraints(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& V,
        const std::vector<EdgeVertexCandidate>& candidates,
        const std::function<bool(double)>& is_active,
        const size_t start_i,
        const size_t end_i,
        Builder& constraint_builder) const;

    void edge_edge_candiates_to_constraints(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& V,
        const std::vector<EdgeEdgeCandidate>& candidates,
        const std::function<bool(double)>& is_active,
        const size_t start_i,
        const size_t end_i,
        Builder& constraint_builder) const;

    void face_vertex_candiates_to_constraints(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& V,
        const std::vector<FaceVertexCandidate>& candidates,
        const std::function<bool(double)>& is_active,
        const size_t start_i,
        const size_t end_i,
        Builder& constraint_builder) const;

    void merge_thread_local_constraints(
        const tbb::enumerable_thread_specific<Builder>& local_storage);
};

} // namespace ipc
