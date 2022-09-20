#pragma once

#include <ipc/collision_mesh.hpp>
#include <ipc/collisions/constraints.hpp>
#include <ipc/friction/friction_constraint.hpp>
#include <ipc/utils/eigen_ext.hpp>

#include <Eigen/Core>
#include <Eigen/Sparse>

namespace ipc {

void construct_friction_constraint_set(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const Constraints& contact_constraint_set,
    double dhat,
    double barrier_stiffness,
    double mu,
    FrictionConstraints& friction_constraint_set);

void construct_friction_constraint_set(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const Constraints& contact_constraint_set,
    double dhat,
    double barrier_stiffness,
    const Eigen::VectorXd& mus,
    FrictionConstraints& friction_constraint_set);

void construct_friction_constraint_set(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V,
    const Constraints& contact_constraint_set,
    double dhat,
    double barrier_stiffness,
    const Eigen::VectorXd& mus,
    const std::function<double(double, double)>& blend_mu,
    FrictionConstraints& friction_constraint_set);

/// @brief Compute the friction potential between two positions.
/// @param[in] mesh The collision mesh.
/// @param[in] V0 Vertex positions at start of time-step (rowwise)
/// @param[in] V1 Current vertex positions (rowwise)
/// @param[in] friction_constraint_set The set of friction constraints.
/// @param[in] epsv_times_h Tolerance for the transition between static and dynamic friction.
template <typename T>
T compute_friction_potential(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V0,
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& V1,
    const FrictionConstraints& friction_constraint_set,
    double epsv_times_h);

/// @brief Compute the gradient of the friction potential wrt V1.
/// @param[in] mesh The collision mesh.
/// @param[in] V0 Vertex positions at start of time-step (rowwise)
/// @param[in] V1 Current vertex positions (rowwise)
/// @param[in] friction_constraint_set The set of friction constraints.
/// @param[in] epsv_times_h Tolerance for the transition between static and dynamic friction.
Eigen::VectorXd compute_friction_potential_gradient(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const FrictionConstraints& friction_constraint_set,
    double epsv_times_h);

/// @brief Compute the Hessian of the friction potential wrt V1.
/// @param[in] mesh The collision mesh.
/// @param[in] V0 Vertex positions at start of time-step (rowwise)
/// @param[in] V1 Current vertex positions (rowwise)
/// @param[in] friction_constraint_set The set of friction constraints.
/// @param[in] epsv_times_h Tolerance for the transition between static and dynamic friction.
Eigen::SparseMatrix<double> compute_friction_potential_hessian(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const FrictionConstraints& friction_constraint_set,
    double epsv_times_h,
    bool project_hessian_to_psd = true);

Eigen::VectorXd compute_friction_force(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& X,
    const Eigen::MatrixXd& Ut,
    const Eigen::MatrixXd& U,
    const FrictionConstraints& friction_constraint_set,
    const double dhat,
    const double barrier_stiffness,
    const double epsv_times_h,
    const double dmin = 0,
    const bool no_mu = false); //< whether to not multiply by mu

inline Eigen::VectorXd compute_friction_force(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& X,
    const Eigen::MatrixXd& U,
    const FrictionConstraints& friction_constraint_set,
    const double dhat,
    const double barrier_stiffness,
    const double epsv_times_h,
    const double dmin = 0,
    const bool no_mu = false) //< whether to not multiply by mu
{
    return compute_friction_force(
        mesh, X, Eigen::MatrixXd::Zero(U.rows(), U.cols()), U,
        friction_constraint_set, dhat, barrier_stiffness, epsv_times_h, dmin,
        no_mu);
}

Eigen::SparseMatrix<double> compute_friction_force_jacobian(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& X,
    const Eigen::MatrixXd& Ut,
    const Eigen::MatrixXd& U,
    const FrictionConstraints& friction_constraint_set,
    const double dhat,
    const double barrier_stiffness,
    const double epsv_times_h,
    const FrictionConstraint::DiffWRT wrt,
    const double dmin = 0);

inline Eigen::SparseMatrix<double> compute_friction_force_jacobian(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& X,
    const Eigen::MatrixXd& U,
    const FrictionConstraints& friction_constraint_set,
    const double dhat,
    const double barrier_stiffness,
    const double epsv_times_h,
    const FrictionConstraint::DiffWRT wrt,
    const double dmin = 0)
{
    return compute_friction_force_jacobian(
        mesh, X, Eigen::MatrixXd::Zero(U.rows(), U.cols()), U,
        friction_constraint_set, dhat, barrier_stiffness, epsv_times_h, wrt,
        dmin);
}

} // namespace ipc

#include "friction.tpp"
