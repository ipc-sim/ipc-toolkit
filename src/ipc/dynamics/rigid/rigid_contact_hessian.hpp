#pragma once

#include <ipc/dynamics/rigid/pose.hpp>
#include <ipc/dynamics/rigid/rigid_bodies.hpp>
#include <ipc/utils/eigen_ext.hpp>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <vector>

namespace ipc::rigid {

/// @brief Assemble a collision potential's Hessian in rigid DOFs, projecting
/// each contact's local block to PSD — a 1-to-1 port of the original Rigid IPC
/// `apply_chain_rule` (Ferguson et al. 2021).
///
/// For each collision the local rigid-DOF Hessian spanning its ≤2 bodies is
///   gn_scale · Jᵢᵀ H_stencil Jᵢ  +  geom_scale · Σ_v (d²Vᵥ/dx²)·grad_stencil,v
/// (Gauss–Newton term with the two-body coupling, plus the indefinite
/// second-order geometric term on the diagonal body blocks). That completed
/// local block is projected with @p project and scattered into the global
/// sparse Hessian.
///
/// This is a 1-to-1 port of Rigid IPC: the interior vertex-space stencil
/// Hessian is left UNPROJECTED, and only the completed rigid-DOF block is
/// projected — *after* adding the indefinite geometric term. Projecting there
/// (rather than the vertex-space block before the geometric term) is what keeps
/// Newton conditioned on stiff, high-mass-ratio contact (e.g. dense piles),
/// cutting iteration counts from hundreds to a handful.
///
/// Defined for (NormalCollisions, BarrierPotential) and
/// (TangentialCollisions, FrictionPotential) via explicit instantiation.
///
/// @param bodies The rigid bodies (collision mesh + per-body transforms).
/// @param poses The current pose of each body.
/// @param collisions The collision set.
/// @param potential The matching potential.
/// @param X Vertex-space quantity the stencil DOFs are read from: world
///          vertices for the barrier, velocities for friction.
/// @param gn_scale Scale applied to the Gauss–Newton term Jᵀ H J.
/// @param geom_scale Scale applied to the geometric term Σ d²V/dx²·grad.
/// @param project PSD projection of each per-contact local block (CLAMP mirrors
///                Rigid IPC; NONE disables projection).
/// @return The potential Hessian in rigid DOFs (size ndof × ndof).
///
/// @note We assemble with plain triplets here instead of going through
/// Potential::assemble_hessian and a HessianAssembler. Those backends scatter
/// per-vertex dim×dim blocks, but our blocks are per body (6×6 in 3D), which
/// the MeshFEMSparse backend does not support. The tradeoff is small: the
/// result is only (6·#bodies)², so assembly is cheap next to the per-contact
/// derivative evaluation, and we never build the large vertex-space Hessian
/// that the block backend was designed to speed up.
template <class TCollisions, class TPotential>
Eigen::SparseMatrix<double> per_contact_rigid_hessian(
    const RigidBodies& bodies,
    const std::vector<Pose>& poses,
    const TCollisions& collisions,
    const TPotential& potential,
    Eigen::ConstRef<Eigen::MatrixXd> X,
    const double gn_scale,
    const double geom_scale,
    const PSDProjectionMethod project);

} // namespace ipc::rigid
