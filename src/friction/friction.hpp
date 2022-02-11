#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <ipc/collision_mesh.hpp>
#include <ipc/collision_constraint.hpp>
#include <ipc/friction/friction_constraint.hpp>
#include <ipc/utils/eigen_ext.hpp>

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

/// @brief Compute the friction potential between to positions.
///
/// @param V0 Vertex positions at start of time-step (rowwise)
/// @param V1 Current vertex positions (rowwise)
/// @param E  Edge vertex indicies
/// @param F  Face vertex indicies (empty in 2D)
/// @param friction_constraint_set
/// @param epsv_times_h
template <typename T>
T compute_friction_potential(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V0,
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& V1,
    const FrictionConstraints& friction_constraint_set,
    double epsv_times_h);

Eigen::VectorXd compute_friction_potential_gradient(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const FrictionConstraints& friction_constraint_set,
    double epsv_times_h);

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
    const double dmin = 0);

inline Eigen::VectorXd compute_friction_force(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& X,
    const Eigen::MatrixXd& U,
    const FrictionConstraints& friction_constraint_set,
    const double dhat,
    const double barrier_stiffness,
    const double epsv_times_h,
    const double dmin = 0)
{
    return compute_friction_force(
        mesh, X, Eigen::MatrixXd(), U, friction_constraint_set, dhat,
        barrier_stiffness, epsv_times_h, dmin);
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
        mesh, X, Eigen::MatrixXd(), U, friction_constraint_set, dhat,
        barrier_stiffness, epsv_times_h, wrt, dmin);
}

} // namespace ipc

#include "friction.tpp"
