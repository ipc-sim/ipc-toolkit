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
        double dhat,
        double barrier_stiffness,
        const Eigen::VectorXd& mus,
        const std::function<double(double, double)>& blend_mu =
            default_blend_mu);

    // ------------------------------------------------------------------------

    /// @brief Compute the friction dissapative potential from the given velocity.
    /// @param mesh The collision mesh.
    /// @param velocity Current vertex velocity (rowwise).
    /// @return The friction dissapative potential.
    template <typename T>
    T compute_potential(
        const CollisionMesh& mesh,
        const MatrixX<T>& velocity,
        double epsv_times_h) const;

    /// @brief Compute the gradient of the friction dissapative potential wrt the velocity.
    /// @param mesh The collision mesh.
    /// @param velocity Current vertex velocity (rowwise).
    /// @return The gradient of the friction dissapative potential wrt the velocity.
    Eigen::VectorXd compute_potential_gradient(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& velocity,
        double epsv_times_h) const;

    /// @brief Compute the Hessian of the friction dissapative potential wrt the velocity.
    /// @param mesh The collision mesh.
    /// @param velocity Current vertex velocity (rowwise).
    /// @param project_hessian_to_psd If true, project the Hessian to be positive semi-definite.
    /// @return The Hessian of the friction dissapative potential wrt the velocity.
    Eigen::SparseMatrix<double> compute_potential_hessian(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& velocity,
        double epsv_times_h,
        bool project_hessian_to_psd = false) const;

    // ------------------------------------------------------------------------

    /// @param no_mu whether to not multiply by mu
    Eigen::VectorXd compute_force(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& X,
        const Eigen::MatrixXd& Ut,
        const Eigen::MatrixXd& U,
        const double dhat,
        const double barrier_stiffness,
        const double epsv_times_h,
        const double dmin = 0,
        const bool no_mu = false) const;

    /// @param no_mu whether to not multiply by mu
    Eigen::VectorXd compute_force(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& X,
        const Eigen::MatrixXd& U,
        const double dhat,
        const double barrier_stiffness,
        const double epsv_times_h,
        const double dmin = 0,
        const bool no_mu = false) const
    {
        return compute_force(
            mesh, X, Eigen::MatrixXd::Zero(U.rows(), U.cols()), U, dhat,
            barrier_stiffness, epsv_times_h, dmin, no_mu);
    }

    Eigen::SparseMatrix<double> compute_force_jacobian(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& X,
        const Eigen::MatrixXd& Ut,
        const Eigen::MatrixXd& U,
        const double dhat,
        const double barrier_stiffness,
        const double epsv_times_h,
        const FrictionConstraint::DiffWRT wrt,
        const double dmin = 0) const;

    Eigen::SparseMatrix<double> compute_force_jacobian(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& X,
        const Eigen::MatrixXd& U,
        const double dhat,
        const double barrier_stiffness,
        const double epsv_times_h,
        const FrictionConstraint::DiffWRT wrt,
        const double dmin = 0) const
    {
        return compute_force_jacobian(
            mesh, X, Eigen::MatrixXd::Zero(U.rows(), U.cols()), U, dhat,
            barrier_stiffness, epsv_times_h, wrt, dmin);
    }

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

#include "friction_constraints.tpp"
