#pragma once

#include <ipc/dynamics/affine/pose.hpp>
#include <ipc/dynamics/rigid/rigid_bodies.hpp>

namespace ipc::affine {

/// @brief Extract the affine pose of body i from the stacked DOF vector.
/// DOF layout per body: [p (dim); vec(A) column-major (dim²)].
affine::Pose
dof_to_pose(Eigen::ConstRef<Eigen::VectorXd> x, const size_t i, const int dim);

/// @brief Extract all affine poses from the stacked DOF vector.
std::vector<affine::Pose>
dof_to_poses(Eigen::ConstRef<Eigen::VectorXd> x, const int dim);

/// @brief Compute the world-space collision mesh vertices from the affine DOFs.
Eigen::MatrixXd
vertices(const rigid::RigidBodies& bodies, Eigen::ConstRef<Eigen::VectorXd> x);

/// @brief The block-diagonal Jacobian dV/dx of the collision mesh vertices
/// with respect to the affine DOFs.
/// @note V(x) is linear in the affine DOFs, so this is constant and the
/// barrier chain rule is simply Jᵀ g / Jᵀ H J with no second-order term.
Eigen::SparseMatrix<double> affine_jacobian(const rigid::RigidBodies& bodies);

} // namespace ipc::affine
