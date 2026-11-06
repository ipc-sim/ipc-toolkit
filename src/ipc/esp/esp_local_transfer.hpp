// Internal implementation of the experimental local quadrature rule.
#pragma once

#include <ipc/esp/esp_collisions.hpp>

namespace ipc::detail {
// c: collisions, m: mesh, X: vertex positions, p: parameters.
double local_transfer_value(
    const ESPCollisions& c,
    const CollisionMesh& m,
    Eigen::ConstRef<Eigen::MatrixXd> X,
    const ESPParameters& p);
Eigen::VectorXd local_transfer_gradient(
    const ESPCollisions& c,
    const CollisionMesh& m,
    Eigen::ConstRef<Eigen::MatrixXd> X,
    const ESPParameters& p);
Eigen::SparseMatrix<double> local_transfer_hessian(
    const ESPCollisions& c,
    const CollisionMesh& m,
    Eigen::ConstRef<Eigen::MatrixXd> X,
    const ESPParameters& p,
    PSDProjectionMethod psd);

} // namespace ipc::detail
