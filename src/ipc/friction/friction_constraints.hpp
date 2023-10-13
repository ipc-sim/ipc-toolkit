#pragma once

#include <ipc/friction/constraints/friction_constraint.hpp>
#include <ipc/friction/constraints/vertex_vertex.hpp>
#include <ipc/friction/constraints/edge_vertex.hpp>
#include <ipc/friction/constraints/edge_edge.hpp>
#include <ipc/friction/constraints/face_vertex.hpp>

#include <ipc/collision_mesh.hpp>
#include <ipc/collisions/collision_constraints.hpp>
#include <ipc/utils/eigen_ext.hpp>

#include <Eigen/Core>
#include <Eigen/Sparse>

namespace ipc {

class FrictionConstraints {
public:
    FrictionConstraints() { }

    void build(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const CollisionConstraints& contact_constraints,
        double dhat,
        double barrier_stiffness,
        double mu)
    {
        this->build(
            mesh, vertices, contact_constraints, dhat, barrier_stiffness,
            Eigen::VectorXd::Constant(vertices.rows(), mu));
    }

    void build(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const CollisionConstraints& contact_constraints,
        const double dhat,
        const double barrier_stiffness,
        const Eigen::VectorXd& mus,
        const std::function<double(double, double)>& blend_mu =
            default_blend_mu);

    // ------------------------------------------------------------------------

    /// @brief Compute the friction dissapative potential from the given velocity.
    /// @param mesh The collision mesh.
    /// @param velocity Current vertex velocity (rowwise).
    /// @param epsv Mollifier parameter \f$\epsilon_v\f$.
    /// @return The friction dissapative potential.
    double compute_potential(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& velocity,
        const double epsv) const;

    /// @brief Compute the gradient of the friction dissapative potential wrt the velocity.
    /// @param mesh The collision mesh.
    /// @param velocity Current vertex velocity (rowwise).
    /// @param epsv Mollifier parameter \f$\epsilon_v\f$.
    /// @return The gradient of the friction dissapative potential wrt the velocity.
    Eigen::VectorXd compute_potential_gradient(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& velocity,
        const double epsv) const;

    /// @brief Compute the Hessian of the friction dissapative potential wrt the velocity.
    /// @param mesh The collision mesh.
    /// @param velocity Current vertex velocity (rowwise).
    /// @param epsv Mollifier parameter \f$\epsilon_v\f$.
    /// @param project_hessian_to_psd If true, project the Hessian to be positive semi-definite.
    /// @return The Hessian of the friction dissapative potential wrt the velocity.
    Eigen::SparseMatrix<double> compute_potential_hessian(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& velocity,
        const double epsv,
        const bool project_hessian_to_psd = false) const;

    // ------------------------------------------------------------------------

    /// @brief Compute the friction force from the given velocity.
    /// @param mesh The collision mesh.
    /// @param rest_positions Rest positions of the vertices (rowwise)
    /// @param lagged_displacements Previous displacements of the vertices (rowwise)
    /// @param velocities Current displacements of the vertices (rowwise)
    /// @param dhat Barrier activation distance.
    /// @param barrier_stiffness Barrier stiffness.
    /// @param epsv Mollifier parameter \f$\epsilon_v\f$.
    /// @param dmin Minimum distance to use for the barrier.
    /// @param no_mu whether to not multiply by mu
    /// @return The friction force.
    Eigen::VectorXd compute_force(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& rest_positions,
        const Eigen::MatrixXd& lagged_displacements,
        const Eigen::MatrixXd& velocities,
        const double dhat,
        const double barrier_stiffness,
        const double epsv,
        const double dmin = 0,
        const bool no_mu = false) const;

    /// @brief Compute the Jacobian of the friction force wrt the velocity.
    /// @param mesh The collision mesh.
    /// @param rest_positions Rest positions of the vertices (rowwise)
    /// @param lagged_displacements Previous displacements of the vertices (rowwise)
    /// @param velocities Current displacements of the vertices (rowwise)
    /// @param dhat Barrier activation distance.
    /// @param barrier_stiffness Barrier stiffness.
    /// @param epsv Mollifier parameter \f$\epsilon_v\f$.
    /// @param wrt The variable to take the derivative with respect to.
    /// @param dmin Minimum distance to use for the barrier.
    /// @return The Jacobian of the friction force wrt the velocity.
    Eigen::SparseMatrix<double> compute_force_jacobian(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& rest_positions,
        const Eigen::MatrixXd& lagged_displacements,
        const Eigen::MatrixXd& velocities,
        const double dhat,
        const double barrier_stiffness,
        const double epsv,
        const FrictionConstraint::DiffWRT wrt,
        const double dmin = 0) const;

    // ------------------------------------------------------------------------

    /// @brief Get the number of friction constraints.
    size_t size() const;

    /// @brief Get if the friction constraints are empty.
    bool empty() const;

    /// @brief Clear the friction constraints.
    void clear();

    /// @brief Get a reference to constriant idx.
    /// @param idx The index of the constraint.
    /// @return A reference to the constraint.
    FrictionConstraint& operator[](const size_t idx);

    /// @brief Get a const reference to constriant idx.
    /// @param idx The index of the constraint.
    /// @return A const reference to the constraint.
    const FrictionConstraint& operator[](const size_t idx) const;

    static double default_blend_mu(double mu0, double mu1)
    {
        // return mu0 * mu1;
        // return std::min(mu0, mu1);
        // return std::max(mu0, mu1);
        return (mu0 + mu1) / 2;
    }

public:
    std::vector<VertexVertexFrictionConstraint> vv_constraints;
    std::vector<EdgeVertexFrictionConstraint> ev_constraints;
    std::vector<EdgeEdgeFrictionConstraint> ee_constraints;
    std::vector<FaceVertexFrictionConstraint> fv_constraints;
};

} // namespace ipc
