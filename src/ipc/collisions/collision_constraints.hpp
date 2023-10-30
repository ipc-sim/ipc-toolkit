#pragma once

#include <ipc/collision_mesh.hpp>
#include <ipc/collisions/collision_constraint.hpp>
#include <ipc/collisions/vertex_vertex.hpp>
#include <ipc/collisions/edge_vertex.hpp>
#include <ipc/collisions/edge_edge.hpp>
#include <ipc/collisions/face_vertex.hpp>
#include <ipc/collisions/plane_vertex.hpp>
#include <ipc/broad_phase/broad_phase.hpp>
#include <ipc/candidates/candidates.hpp>

#include <Eigen/Core>

#include <vector>

namespace ipc {

class CollisionConstraints {
public:
    CollisionConstraints() { }

    /// @brief Initialize the set of constraints used to compute the barrier potential.
    /// @param mesh The collision mesh.
    /// @param vertices Vertices of the collision mesh.
    /// @param dhat The activation distance of the barrier.
    /// @param dmin Minimum distance.
    /// @param broad_phase_method Broad-phase method to use.
    void build(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const double dhat,
        const double dmin = 0,
        const BroadPhaseMethod broad_phase_method = DEFAULT_BROAD_PHASE_METHOD);

    /// @brief Initialize the set of constraints used to compute the barrier potential.
    /// @param candidates Distance candidates from which the constraint set is built.
    /// @param mesh The collision mesh.
    /// @param vertices Vertices of the collision mesh.
    /// @param dhat The activation distance of the barrier.
    /// @param  dmin  Minimum distance.
    void build(
        const Candidates& candidates,
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const double dhat,
        const double dmin = 0);

    // ------------------------------------------------------------------------

    /// @brief Compute the barrier potential for a given constraint set.
    /// @param mesh The collision mesh.
    /// @param vertices Vertices of the collision mesh.
    /// @param dhat The activation distance of the barrier.
    /// @returns The sum of all barrier potentials (not scaled by the barrier stiffness).
    double compute_potential(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const double dhat) const;

    /// @brief Compute the gradient of the barrier potential.
    /// @param mesh The collision mesh.
    /// @param vertices Vertices of the collision mesh.
    /// @param dhat The activation distance of the barrier.
    /// @returns The gradient of all barrier potentials (not scaled by the barrier stiffness). This will have a size of |vertices|.
    Eigen::VectorXd compute_potential_gradient(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const double dhat) const;

    /// @brief Compute the hessian of the barrier potential.
    /// @param mesh The collision mesh.
    /// @param vertices Vertices of the collision mesh.
    /// @param dhat The activation distance of the barrier.
    /// @param project_hessian_to_psd Make sure the hessian is positive semi-definite.
    /// @returns The hessian of all barrier potentials (not scaled by the barrier stiffness). This will have a size of |vertices|x|vertices|.
    Eigen::SparseMatrix<double> compute_potential_hessian(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const double dhat,
        const bool project_hessian_to_psd = false) const;

    // ------------------------------------------------------------------------

    /// @brief Compute the barrier shape derivative.
    /// @param mesh The collision mesh.
    /// @param vertices Vertices of the collision mesh.
    /// @param dhat The activation distance of the barrier.
    /// @throws std::runtime_error If the collision constraints were not built with shape derivatives enabled.
    /// @returns The derivative of the force with respect to X, the rest vertices.
    Eigen::SparseMatrix<double> compute_shape_derivative(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const double dhat) const;

    /// @brief Computes the minimum distance between any non-adjacent elements.
    /// @param mesh The collision mesh.
    /// @param vertices Vertices of the collision mesh.
    /// @returns The minimum distance between any non-adjacent elements.
    double compute_minimum_distance(
        const CollisionMesh& mesh, const Eigen::MatrixXd& vertices) const;

    // ------------------------------------------------------------------------

    /// @brief Get the number of collision constraints.
    size_t size() const;

    /// @brief Get if the collision constraints are empty.
    bool empty() const;

    /// @brief Clear the collision constraints.
    void clear();

    /// @brief Get a reference to constriant idx.
    /// @param idx The index of the constraint.
    /// @return A reference to the constraint.
    CollisionConstraint& operator[](size_t idx);

    /// @brief Get a const reference to constriant idx.
    /// @param idx The index of the constraint.
    /// @return A const reference to the constraint.
    const CollisionConstraint& operator[](size_t idx) const;

    /// @brief Get if the constraint at idx is a vertex-vertex constraint.
    /// @param idx The index of the constraint.
    /// @return If the constraint at idx is a vertex-vertex constraint.
    bool is_vertex_vertex(size_t idx) const;

    /// @brief Get if the constraint at idx is an edge-vertex constraint.
    /// @param idx The index of the constraint.
    /// @return If the constraint at idx is an edge-vertex constraint.
    bool is_edge_vertex(size_t idx) const;

    /// @brief Get if the constraint at idx is an edge-edge constraint.
    /// @param idx The index of the constraint.
    /// @return If the constraint at idx is an edge-edge constraint.
    bool is_edge_edge(size_t idx) const;

    /// @brief Get if the constraint at idx is an face-vertex constraint.
    /// @param idx The index of the constraint.
    /// @return If the constraint at idx is an face-vertex constraint.
    bool is_face_vertex(size_t idx) const;

    /// @brief Get if the constraint at idx is an plane-vertex constraint.
    /// @param idx The index of the constraint.
    /// @return If the constraint at idx is an plane-vertex constraint.
    bool is_plane_vertex(size_t idx) const;

    /// @brief Get if the collision constraints should use the convergent formulation.
    /// @note If not empty, this is the current value not necessarily the value used to build the constraints.
    /// @return If the collision constraints should use the convergent formulation.
    bool use_convergent_formulation() const
    {
        return m_use_convergent_formulation;
    }

    /// @brief Set if the collision constraints should use the convergent formulation.
    /// @warning This must be set before the constraints are built.
    /// @param use_convergent_formulation If the collision constraints should use the convergent formulation.
    void set_use_convergent_formulation(const bool use_convergent_formulation);

    /// @brief Get if the collision constraints are using the convergent formulation.
    /// @note If not empty, this is the current value not necessarily the value used to build the constraints.
    /// @return If the collision constraints are using the convergent formulation.
    bool are_shape_derivatives_enabled() const
    {
        return m_are_shape_derivatives_enabled;
    }

    /// @brief Set if the collision constraints should enable shape derivative computation.
    /// @warning This must be set before the constraints are built.
    /// @param are_shape_derivatives_enabled If the collision constraints should enable shape derivative computation.
    void
    set_are_shape_derivatives_enabled(const bool are_shape_derivatives_enabled);

    std::string
    to_string(const CollisionMesh& mesh, const Eigen::MatrixXd& vertices) const;

public:
    std::vector<VertexVertexConstraint> vv_constraints;
    std::vector<EdgeVertexConstraint> ev_constraints;
    std::vector<EdgeEdgeConstraint> ee_constraints;
    std::vector<FaceVertexConstraint> fv_constraints;
    std::vector<PlaneVertexConstraint> pv_constraints;

protected:
    bool m_use_convergent_formulation = false;
    bool m_are_shape_derivatives_enabled = false;
};

} // namespace ipc
