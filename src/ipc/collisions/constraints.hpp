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

#include <vector>

namespace ipc {

class CollisionConstraints {
public:
    CollisionConstraints() { }

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

    // ------------------------------------------------------------------------

    /// @brief Compute the barrier potential for a given constraint set.
    /// @param[in] mesh The collision mesh.
    /// @param[in] V Vertices of the collision mesh.
    /// @param[in] constraint_set The set of constraints.
    /// @param[in] dhat The activation distance of the barrier.
    /// @returns The sum of all barrier potentials (not scaled by the barrier stiffness).
    double compute_potential(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& V,
        const double dhat) const;

    /// @brief Compute the gradient of the barrier potential.
    /// @param[in] mesh The collision mesh.
    /// @param[in] V Vertices of the collision mesh.
    /// @param[in] constraint_set The set of constraints.
    /// @param[in] dhat The activation distance of the barrier.
    /// @returns The gradient of all barrier potentials (not scaled by the barrier stiffness). This will have a size of |V|.
    Eigen::VectorXd compute_potential_gradient(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& V,
        const double dhat) const;

    /// @brief Compute the hessian of the barrier potential.
    /// @param[in] mesh The collision mesh.
    /// @param[in] V Vertices of the collision mesh.
    /// @param[in] constraint_set The set of constraints.
    /// @param[in] dhat The activation distance of the barrier.
    /// @param[in] project_hessian_to_psd Make sure the hessian is positive semi-definite.
    /// @returns The hessian of all barrier potentials (not scaled by the barrier stiffness). This will have a size of |V|x|V|.
    Eigen::SparseMatrix<double> compute_potential_hessian(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& V,
        const double dhat,
        const bool project_hessian_to_psd = true) const;

    // ------------------------------------------------------------------------

    /// @brief Compute the barrier shape derivative.
    /// @param[in] mesh The collision mesh.
    /// @param[in] V Vertices of the collision mesh.
    /// @param[in] constraint_set The set of constraints.
    /// @param[in] dhat The activation distance of the barrier.
    /// @returns The derivative of the force with respect to X, the rest positions.
    Eigen::SparseMatrix<double> compute_shape_derivative(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& V,
        const double dhat) const;

    /// @brief Computes the minimum distance between any non-adjacent elements.
    /// @param[in] mesh The collision mesh.
    /// @param[in] V Vertices of the collision mesh.
    /// @returns The minimum distance between any non-adjacent elements.
    double compute_minimum_distance(
        const CollisionMesh& mesh, const Eigen::MatrixXd& V) const;

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

public:
    std::vector<VertexVertexConstraint> vv_constraints;
    std::vector<EdgeVertexConstraint> ev_constraints;
    std::vector<EdgeEdgeConstraint> ee_constraints;
    std::vector<FaceVertexConstraint> fv_constraints;
    std::vector<PlaneVertexConstraint> pv_constraints;
    bool use_convergent_formulation = false;
    bool compute_shape_derivatives = false;
};

} // namespace ipc
