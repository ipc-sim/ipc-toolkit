#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <ipc/collision_constraint.hpp>
#include <ipc/friction/friction_constraint.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

void construct_friction_constraint_set(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const Constraints& contact_constraint_set,
    double dhat,
    double barrier_stiffness,
    double mu,
    FrictionConstraints& friction_constraint_set);

/// @brief Compute the friction potential between to positions.
///
/// @param V0 Vertex positions at start of time-step (or newton-solve) (rowwise)
/// @param V1 Current vertex positions (rowwise)
/// @param E  Edge vertex indicies
/// @param F  Face vertex indicies (empty in 2D)
/// @param friction_constraint_set
/// @param closest_points
/// @param tangent_bases
/// @param normal_force_magnitudes
/// @param epsv_times_h
/// @param mu Coefficient of friction for all contacts
template <typename T>
T compute_friction_potential(
    const Eigen::MatrixXd& V0,
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& V1,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const FrictionConstraints& friction_constraint_set,
    double epsv_times_h);

Eigen::VectorXd compute_friction_potential_gradient(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const FrictionConstraints& friction_constraint_set,
    double epsv_times_h);

Eigen::SparseMatrix<double> compute_friction_potential_hessian(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const FrictionConstraints& friction_constraint_set,
    double epsv_times_h,
    bool project_hessian_to_psd = true);

///////////////////////////////////////////////////////////////////////////////
// Compute potential and derivatives for a single constraint.

template <
    typename DerivedDP0,
    typename DerivedDP1,
    typename T = typename DerivedDP0::Scalar>
inline T compute_friction_potential(
    const Eigen::MatrixBase<DerivedDP0>& dp0,
    const Eigen::MatrixBase<DerivedDP1>& dp1,
    const VertexVertexFrictionConstraint& friction_constraint,
    double epsv_times_h);

template <
    typename DerivedDP,
    typename DerivedDE0,
    typename DerivedDE1,
    typename T = typename DerivedDP::Scalar>
inline T compute_friction_potential(
    const Eigen::MatrixBase<DerivedDP>& dp,
    const Eigen::MatrixBase<DerivedDE0>& de0,
    const Eigen::MatrixBase<DerivedDE1>& de1,
    const EdgeVertexFrictionConstraint& friction_constraint,
    double epsv_times_h);

template <
    typename DerivedDEA0,
    typename DerivedDEA1,
    typename DerivedDEB0,
    typename DerivedDEB1,
    typename T = typename DerivedDEA0::Scalar>
inline T compute_friction_potential(
    const Eigen::MatrixBase<DerivedDEA0>& dea0,
    const Eigen::MatrixBase<DerivedDEA1>& dea1,
    const Eigen::MatrixBase<DerivedDEB0>& deb0,
    const Eigen::MatrixBase<DerivedDEB1>& deb1,
    const EdgeEdgeFrictionConstraint& friction_constraint,
    double epsv_times_h);

template <
    typename DerivedDP,
    typename DerivedDT0,
    typename DerivedDT1,
    typename DerivedDT2,
    typename T = typename DerivedDP::Scalar>
inline T compute_friction_potential(
    const Eigen::MatrixBase<DerivedDP>& dp,
    const Eigen::MatrixBase<DerivedDT0>& dt0,
    const Eigen::MatrixBase<DerivedDT1>& dt1,
    const Eigen::MatrixBase<DerivedDT2>& dt2,
    const FaceVertexFrictionConstraint& friction_constraint,
    double epsv_times_h);

///////////////////////////////////////////////////////////////////////////////

template <typename DerivedDP0, typename DerivedDP1>
inline VectorMax6d compute_friction_potential_gradient(
    const Eigen::MatrixBase<DerivedDP0>& dp0,
    const Eigen::MatrixBase<DerivedDP1>& dp1,
    const VertexVertexFrictionConstraint& friction_constraint,
    double epsv_times_h);

template <typename DerivedDP, typename DerivedDE0, typename DerivedDE1>
inline VectorMax9d compute_friction_potential_gradient(
    const Eigen::MatrixBase<DerivedDP>& dp,
    const Eigen::MatrixBase<DerivedDE0>& de0,
    const Eigen::MatrixBase<DerivedDE1>& de1,
    const EdgeVertexFrictionConstraint& friction_constraint,
    double epsv_times_h);

template <
    typename DerivedDEA0,
    typename DerivedDEA1,
    typename DerivedDEB0,
    typename DerivedDEB1>
inline VectorMax12d compute_friction_potential_gradient(
    const Eigen::MatrixBase<DerivedDEA0>& dea0,
    const Eigen::MatrixBase<DerivedDEA1>& dea1,
    const Eigen::MatrixBase<DerivedDEB0>& deb0,
    const Eigen::MatrixBase<DerivedDEB1>& deb1,
    const EdgeEdgeFrictionConstraint& friction_constraint,
    double epsv_times_h);

template <
    typename DerivedDP,
    typename DerivedDT0,
    typename DerivedDT1,
    typename DerivedDT2>
inline VectorMax12d compute_friction_potential_gradient(
    const Eigen::MatrixBase<DerivedDP>& dp,
    const Eigen::MatrixBase<DerivedDT0>& dt0,
    const Eigen::MatrixBase<DerivedDT1>& dt1,
    const Eigen::MatrixBase<DerivedDT2>& dt2,
    const FaceVertexFrictionConstraint& friction_constraint,
    double epsv_times_h);

///////////////////////////////////////////////////////////////////////////////

template <typename DerivedDP0, typename DerivedDP1>
inline MatrixMax12d compute_friction_potential_hessian(
    const Eigen::MatrixBase<DerivedDP0>& dp0,
    const Eigen::MatrixBase<DerivedDP1>& dp1,
    const VertexVertexFrictionConstraint& friction_constraint,
    double epsv_times_h,
    bool project_hessian_to_psd = true);

template <typename DerivedDP, typename DerivedDE0, typename DerivedDE1>
inline MatrixMax12d compute_friction_potential_hessian(
    const Eigen::MatrixBase<DerivedDP>& dp,
    const Eigen::MatrixBase<DerivedDE0>& de0,
    const Eigen::MatrixBase<DerivedDE1>& de1,
    const EdgeVertexFrictionConstraint& friction_constraint,
    double epsv_times_h,
    bool project_hessian_to_psd = true);

template <
    typename DerivedDEA0,
    typename DerivedDEA1,
    typename DerivedDEB0,
    typename DerivedDEB1>
inline MatrixMax12d compute_friction_potential_hessian(
    const Eigen::MatrixBase<DerivedDEA0>& dea0,
    const Eigen::MatrixBase<DerivedDEA1>& dea1,
    const Eigen::MatrixBase<DerivedDEB0>& deb0,
    const Eigen::MatrixBase<DerivedDEB1>& deb1,
    const EdgeEdgeFrictionConstraint& friction_constraint,
    double epsv_times_h,
    bool project_hessian_to_psd = true);

template <
    typename DerivedDP,
    typename DerivedDT0,
    typename DerivedDT1,
    typename DerivedDT2>
inline MatrixMax12d compute_friction_potential_hessian(
    const Eigen::MatrixBase<DerivedDP>& dp,
    const Eigen::MatrixBase<DerivedDT0>& dt0,
    const Eigen::MatrixBase<DerivedDT1>& dt1,
    const Eigen::MatrixBase<DerivedDT2>& dt2,
    const FaceVertexFrictionConstraint& friction_constraint,
    double epsv_times_h,
    bool project_hessian_to_psd = true);

} // namespace ipc

#include "friction.tpp"
