#include "collision_stencil.hpp"

namespace ipc {

double CollisionStencil::compute_distance(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces) const
{
    return compute_distance(dof(vertices, edges, faces));
}

VectorMax12d CollisionStencil::compute_distance_gradient(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces) const
{
    return compute_distance_gradient(dof(vertices, edges, faces));
}

MatrixMax12d CollisionStencil::compute_distance_hessian(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces) const
{
    return compute_distance_hessian(dof(vertices, edges, faces));
}

} // namespace ipc