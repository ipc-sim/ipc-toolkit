#pragma once

#include <ipc/collision_mesh.hpp>
#include <ipc/collisions/constraints.hpp>
#include <ipc/friction/friction_constraint.hpp>
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
        const Constraints& contact_constraints,
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
        const Constraints& contact_constraints,
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
        const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& velocity,
        double epsv_times_h);

    /// @brief Compute the gradient of the friction dissapative potential wrt the velocity.
    /// @param mesh The collision mesh.
    /// @param velocity Current vertex velocity (rowwise).
    /// @return The gradient of the friction dissapative potential wrt the velocity.
    Eigen::VectorXd compute_potential_gradient(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& velocity,
        double epsv_times_h);

    /// @brief Compute the Hessian of the friction dissapative potential wrt the velocity.
    /// @param mesh The collision mesh.
    /// @param velocity Current vertex velocity (rowwise).
    /// @param project_hessian_to_psd If true, project the Hessian to be positive semi-definite.
    /// @return The Hessian of the friction dissapative potential wrt the velocity.
    Eigen::SparseMatrix<double> compute_potential_hessian(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& velocity,
        double epsv_times_h,
        bool project_hessian_to_psd = true);

    // ------------------------------------------------------------------------

    Eigen::VectorXd compute_force(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& X,
        const Eigen::MatrixXd& Ut,
        const Eigen::MatrixXd& U,
        const double dhat,
        const double barrier_stiffness,
        const double epsv_times_h,
        const double dmin = 0,
        const bool no_mu = false); //< whether to not multiply by mu

    Eigen::VectorXd compute_force(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& X,
        const Eigen::MatrixXd& U,
        const double dhat,
        const double barrier_stiffness,
        const double epsv_times_h,
        const double dmin = 0,
        const bool no_mu = false) //< whether to not multiply by mu
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
        const double dmin = 0);

    Eigen::SparseMatrix<double> compute_force_jacobian(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& X,
        const Eigen::MatrixXd& U,
        const double dhat,
        const double barrier_stiffness,
        const double epsv_times_h,
        const FrictionConstraint::DiffWRT wrt,
        const double dmin = 0)
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

public:
    std::vector<VertexVertexFrictionConstraint> vv_constraints;
    std::vector<EdgeVertexFrictionConstraint> ev_constraints;
    std::vector<EdgeEdgeFrictionConstraint> ee_constraints;
    std::vector<FaceVertexFrictionConstraint> fv_constraints;

private:
    static double default_blend_mu(double mu0, double mu1)
    {
        // return mu0 * mu1;
        // return std::min(mu0, mu1);
        // return std::max(mu0, mu1);
        return (mu0 + mu1) / 2;
    }
};

// void construct_friction_constraints(
//     const CollisionMesh& mesh,
//     const Eigen::MatrixXd& V,
//     const Constraints& contact_constraints,
//     double dhat,
//     double barrier_stiffness,
//     double mu,
//     FrictionConstraints& friction_constraints);

// void construct_friction_constraints(
//     const CollisionMesh& mesh,
//     const Eigen::MatrixXd& V,
//     const Constraints& contact_constraints,
//     double dhat,
//     double barrier_stiffness,
//     const Eigen::VectorXd& mus,
//     FrictionConstraints& friction_constraints);

// void construct_friction_constraints(
//     const CollisionMesh& mesh,
//     const Eigen::MatrixXd& V,
//     const Constraints& contact_constraints,
//     double dhat,
//     double barrier_stiffness,
//     const Eigen::VectorXd& mus,
//     const std::function<double(double, double)>& blend_mu,
//     FrictionConstraints& friction_constraints);

// /// @brief Compute the friction potential between two positions.
// /// @param[in] mesh The collision mesh.
// /// @param[in] V0 Vertex positions at start of time-step (rowwise)
// /// @param[in] V1 Current vertex positions (rowwise)
// /// @param[in] friction_constraints The set of friction constraints.
// /// @param[in] epsv_times_h Tolerance for the transition between static and
// dynamic friction. template <typename T> T compute_friction_potential(
//     const CollisionMesh& mesh,
//     const Eigen::MatrixXd& V0,
//     const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& V1,
//     const FrictionConstraints& friction_constraints,
//     double epsv_times_h);

// /// @brief Compute the gradient of the friction potential wrt V1.
// /// @param[in] mesh The collision mesh.
// /// @param[in] V0 Vertex positions at start of time-step (rowwise)
// /// @param[in] V1 Current vertex positions (rowwise)
// /// @param[in] friction_constraints The set of friction constraints.
// /// @param[in] epsv_times_h Tolerance for the transition between static and
// dynamic friction. Eigen::VectorXd compute_friction_potential_gradient(
//     const CollisionMesh& mesh,
//     const Eigen::MatrixXd& V0,
//     const Eigen::MatrixXd& V1,
//     const FrictionConstraints& friction_constraints,
//     double epsv_times_h);

// /// @brief Compute the Hessian of the friction potential wrt V1.
// /// @param[in] mesh The collision mesh.
// /// @param[in] V0 Vertex positions at start of time-step (rowwise)
// /// @param[in] V1 Current vertex positions (rowwise)
// /// @param[in] friction_constraints The set of friction constraints.
// /// @param[in] epsv_times_h Tolerance for the transition between static and
// dynamic friction. Eigen::SparseMatrix<double>
// compute_friction_potential_hessian(
//     const CollisionMesh& mesh,
//     const Eigen::MatrixXd& V0,
//     const Eigen::MatrixXd& V1,
//     const FrictionConstraints& friction_constraints,
//     double epsv_times_h,
//     bool project_hessian_to_psd = true);

// Eigen::VectorXd compute_friction_force(
//     const CollisionMesh& mesh,
//     const Eigen::MatrixXd& X,
//     const Eigen::MatrixXd& Ut,
//     const Eigen::MatrixXd& U,
//     const FrictionConstraints& friction_constraints,
//     const double dhat,
//     const double barrier_stiffness,
//     const double epsv_times_h,
//     const double dmin = 0,
//     const bool no_mu = false); //< whether to not multiply by mu

// inline Eigen::VectorXd compute_friction_force(
//     const CollisionMesh& mesh,
//     const Eigen::MatrixXd& X,
//     const Eigen::MatrixXd& U,
//     const FrictionConstraints& friction_constraints,
//     const double dhat,
//     const double barrier_stiffness,
//     const double epsv_times_h,
//     const double dmin = 0,
//     const bool no_mu = false) //< whether to not multiply by mu
// {
//     return compute_friction_force(
//         mesh, X, Eigen::MatrixXd::Zero(U.rows(), U.cols()), U,
//         friction_constraints, dhat, barrier_stiffness, epsv_times_h, dmin,
//         no_mu);
// }

// Eigen::SparseMatrix<double> compute_friction_force_jacobian(
//     const CollisionMesh& mesh,
//     const Eigen::MatrixXd& X,
//     const Eigen::MatrixXd& Ut,
//     const Eigen::MatrixXd& U,
//     const FrictionConstraints& friction_constraints,
//     const double dhat,
//     const double barrier_stiffness,
//     const double epsv_times_h,
//     const FrictionConstraint::DiffWRT wrt,
//     const double dmin = 0);

// inline Eigen::SparseMatrix<double> compute_friction_force_jacobian(
//     const CollisionMesh& mesh,
//     const Eigen::MatrixXd& X,
//     const Eigen::MatrixXd& U,
//     const FrictionConstraints& friction_constraints,
//     const double dhat,
//     const double barrier_stiffness,
//     const double epsv_times_h,
//     const FrictionConstraint::DiffWRT wrt,
//     const double dmin = 0)
// {
//     return compute_friction_force_jacobian(
//         mesh, X, Eigen::MatrixXd::Zero(U.rows(), U.cols()), U,
//         friction_constraints, dhat, barrier_stiffness, epsv_times_h, wrt,
//         dmin);
// }

} // namespace ipc

#include "friction.tpp"
