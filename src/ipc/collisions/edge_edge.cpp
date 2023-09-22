#include "edge_edge.hpp"

#include <ipc/barrier/barrier.hpp>
#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/edge_edge_mollifier.hpp>

namespace ipc {

EdgeEdgeConstraint::EdgeEdgeConstraint(
    long edge0_id,
    long edge1_id,
    double eps_x,
    const EdgeEdgeDistanceType dtype)
    : EdgeEdgeCandidate(edge0_id, edge1_id)
    , eps_x(eps_x)
    , dtype(dtype)
{
}

EdgeEdgeConstraint::EdgeEdgeConstraint(
    const EdgeEdgeCandidate& candidate,
    double eps_x,
    const EdgeEdgeDistanceType dtype)
    : EdgeEdgeCandidate(candidate)
    , eps_x(eps_x)
    , dtype(dtype)
{
}

double EdgeEdgeConstraint::compute_potential(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const double dhat) const
{
    return edge_edge_mollifier(
               vertices.row(edges(edge0_id, 0)),
               vertices.row(edges(edge0_id, 1)),
               vertices.row(edges(edge1_id, 0)),
               vertices.row(edges(edge1_id, 1)), eps_x)
        * CollisionConstraint::compute_potential(vertices, edges, faces, dhat);
}

VectorMax12d EdgeEdgeConstraint::compute_potential_gradient(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const double dhat) const
{
    const auto& [ea0, ea1, eb0, eb1] = this->vertices(vertices, edges, faces);

    // b(d(x))
    const double barrier =
        CollisionConstraint::compute_potential(vertices, edges, faces, dhat);
    // ∇ b(d(x))
    const VectorMax12d barrier_grad =
        CollisionConstraint::compute_potential_gradient(
            vertices, edges, faces, dhat);

    // m(x)
    const double mollifier = edge_edge_mollifier(ea0, ea1, eb0, eb1, eps_x);
    // ∇ m(x)
    const Vector12d mollifier_grad =
        edge_edge_mollifier_gradient(ea0, ea1, eb0, eb1, eps_x);

    // ∇[m(x) * b(d(x))] = ∇m(x)) * b(d(x)) + m(x) * ∇ b(d(x))
    return mollifier_grad * barrier + mollifier * barrier_grad;
}

MatrixMax12d EdgeEdgeConstraint::compute_potential_hessian(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const double dhat,
    const bool project_hessian_to_psd) const
{
    const auto& [ea0, ea1, eb0, eb1] = this->vertices(vertices, edges, faces);

    // b(d(x))
    const double barrier =
        CollisionConstraint::compute_potential(vertices, edges, faces, dhat);
    // ∇ b(d(x))
    const Vector12d barrier_grad =
        CollisionConstraint::compute_potential_gradient(
            vertices, edges, faces, dhat);
    // ∇² b(d(x))
    const Matrix12d barrier_hess =
        CollisionConstraint::compute_potential_hessian(
            vertices, edges, faces, dhat, /*project_hessian_to_psd=*/false);

    // m(x)
    const double mollifier = edge_edge_mollifier(ea0, ea1, eb0, eb1, eps_x);
    // ∇ m(x)
    const Vector12d mollifier_grad =
        edge_edge_mollifier_gradient(ea0, ea1, eb0, eb1, eps_x);
    // ∇² m(x)
    const Matrix12d mollifier_hess =
        edge_edge_mollifier_hessian(ea0, ea1, eb0, eb1, eps_x);

    // ∇²[m(x) * b(d(x))] = ∇[∇m(x) * b(d(x)) + m(x) * ∇b(d(x))]
    //                    = ∇²m(x) * b(d(x)) + ∇b(d(x)) * ∇m(x)ᵀ
    //                      + ∇m(x) * ∇b(d(x))ᵀ + m(x) * ∇²b(d(x))
    const Matrix12d grad_b_grad_m = barrier_grad * mollifier_grad.transpose();

    const Matrix12d hess = mollifier_hess * barrier + grad_b_grad_m
        + grad_b_grad_m.transpose() + mollifier * barrier_hess;

    return project_hessian_to_psd ? project_to_psd(hess) : hess;
}

} // namespace ipc
