// Internal implementation of the experimental local quadrature rule.
#pragma once

#include <ipc/esp/esp_collisions.hpp>

namespace ipc::detail {
double local_transfer_value(
    const ESPCollisions&,
    const CollisionMesh&,
    Eigen::ConstRef<Eigen::MatrixXd>,
    const ESPParameters&);
Eigen::VectorXd local_transfer_gradient(
    const ESPCollisions&,
    const CollisionMesh&,
    Eigen::ConstRef<Eigen::MatrixXd>,
    const ESPParameters&);
Eigen::SparseMatrix<double> local_transfer_hessian(
    const ESPCollisions&,
    const CollisionMesh&,
    Eigen::ConstRef<Eigen::MatrixXd>,
    const ESPParameters&,
    PSDProjectionMethod);

// Lagged normalization: the transferred terms only, with per-source-edge
// coefficients sum_{e ∋ ea} w_e c_e s_e / D_e (see esp_potential.cpp).
double lagged_transfer_value(
    const ESPCollisions&,
    const CollisionMesh&,
    Eigen::ConstRef<Eigen::MatrixXd>,
    const ESPParameters&,
    const Eigen::VectorXd& edge_coef);
Eigen::VectorXd lagged_transfer_gradient(
    const ESPCollisions&,
    const CollisionMesh&,
    Eigen::ConstRef<Eigen::MatrixXd>,
    const ESPParameters&,
    const Eigen::VectorXd& edge_coef);
std::vector<Eigen::Triplet<double>> lagged_transfer_hessian(
    const ESPCollisions&,
    const CollisionMesh&,
    Eigen::ConstRef<Eigen::MatrixXd>,
    const ESPParameters&,
    const Eigen::VectorXd& edge_coef,
    PSDProjectionMethod);
} // namespace ipc::detail
