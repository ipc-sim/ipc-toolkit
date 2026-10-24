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
} // namespace ipc::detail
