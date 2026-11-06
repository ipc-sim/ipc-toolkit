#pragma once

#include <ipc/dynamics/rigid/rigid_bodies.hpp>

namespace ipc::affine {

class OrthogonalityPotential {
public:
    OrthogonalityPotential(const double stiffness) : stiffness(stiffness) { }

    // ---- Cumulative functions -----------------------------------------------

    double operator()(
        const rigid::RigidBodies& bodies,
        Eigen::ConstRef<Eigen::VectorXd> x) const;

    Eigen::VectorXd gradient(
        const rigid::RigidBodies& bodies,
        Eigen::ConstRef<Eigen::VectorXd> x) const;

    Eigen::SparseMatrix<double> hessian(
        const rigid::RigidBodies& bodies,
        Eigen::ConstRef<Eigen::VectorXd> x,
        const PSDProjectionMethod project_hessian_to_psd =
            PSDProjectionMethod::NONE) const;

    // ---- Per-body functions -------------------------------------------------

    double operator()(
        const rigid::RigidBody& body, Eigen::ConstRef<VectorMax12d> x) const;

    VectorMax12d gradient(
        const rigid::RigidBody& body, Eigen::ConstRef<VectorMax12d> x) const;

    MatrixMax12d hessian(
        const rigid::RigidBody& body,
        Eigen::ConstRef<VectorMax12d> x,
        const PSDProjectionMethod project_hessian_to_psd =
            PSDProjectionMethod::NONE) const;

    double stiffness;
};

} // namespace ipc::affine