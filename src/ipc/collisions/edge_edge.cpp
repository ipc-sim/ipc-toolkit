#include "edge_edge.hpp"

#include <ipc/barrier/barrier.hpp>
#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/edge_edge_mollifier.hpp>

namespace ipc {

EdgeEdgeConstraint::EdgeEdgeConstraint(
    const long _edge0_id,
    const long _edge1_id,
    const double _eps_x,
    const EdgeEdgeDistanceType _dtype)
    : EdgeEdgeCandidate(_edge0_id, _edge1_id)
    , eps_x(_eps_x)
    , dtype(_dtype)
{
}

EdgeEdgeConstraint::EdgeEdgeConstraint(
    const EdgeEdgeCandidate& candidate,
    const double _eps_x,
    const EdgeEdgeDistanceType _dtype)
    : EdgeEdgeCandidate(candidate)
    , eps_x(_eps_x)
    , dtype(_dtype)
{
}

EdgeEdgeConstraint::EdgeEdgeConstraint(
    const long _edge0_id,
    const long _edge1_id,
    const double _eps_x,
    const double _weight,
    const Eigen::SparseVector<double>& _weight_gradient,
    const EdgeEdgeDistanceType _dtype)
    : EdgeEdgeCandidate(_edge0_id, _edge1_id)
    , CollisionConstraint(_weight, _weight_gradient)
    , eps_x(_eps_x)
    , dtype(_dtype)
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

MatrixMax12d EdgeEdgeConstraint::compute_shape_derivative_second_term(
    const Eigen::MatrixXd& rest_positions, // = x
    const Eigen::MatrixXd& vertices,       // = x + u
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const double dhat) const
{
    const auto& [ea0_rest, ea1_rest, eb0_rest, eb1_rest] =
        this->vertices(rest_positions, edges, faces);
    const auto& [ea0, ea1, eb0, eb1] = this->vertices(vertices, edges, faces);

    // b(d(x+u))
    const double barrier =
        CollisionConstraint::compute_potential(vertices, edges, faces, dhat);
    // ∇ᵤ b(d(x+u))
    const Vector12d barrier_grad =
        CollisionConstraint::compute_potential_gradient(
            vertices, edges, faces, dhat);
    // ∇ᵤ² b(d(x+u))
    const Matrix12d barrier_hess =
        CollisionConstraint::compute_potential_hessian(
            vertices, edges, faces, dhat, /*project_hessian_to_psd=*/false);

    // ε(x)
    const double _eps_x =
        edge_edge_mollifier_threshold(ea1_rest, ea0_rest, eb0_rest, eb1_rest);

    // m(x,u)
    const double mollifier = edge_edge_mollifier(ea0, ea1, eb0, eb1, _eps_x);
    // ∇ᵤ m(x,u)
    const Vector12d mollifier_gradu =
        edge_edge_mollifier_gradient(ea0, ea1, eb0, eb1, _eps_x);
    // ∇ₓ m(x,u)
    const Vector12d mollifier_gradx = edge_edge_mollifier_gradient_wrt_x(
        ea0_rest, ea1_rest, eb0_rest, eb1_rest, ea0, ea1, eb0, eb1);
    // ∇ₓ∇ᵤ m(x,u)
    const Matrix12d mollifier_jacobian =
        edge_edge_mollifier_gradient_jacobian_wrt_x(
            ea0_rest, ea1_rest, eb0_rest, eb1_rest, ea0, ea1, eb0, eb1);

    // Only compute the second term of the shape derivative
    // ∇ₓ (b ∇ᵤm + m ∇ᵤb) = b ∇ₓ∇ᵤm + (∇ₓm)(∇ᵤb)ᵀ + (∇ᵤb)(∇ᵤm)ᵀ + m ∇ᵤ²b
    return barrier * mollifier_jacobian
        + mollifier_gradx * barrier_grad.transpose()
        + barrier_grad * mollifier_gradu.transpose() + mollifier * barrier_hess;
}

bool EdgeEdgeConstraint::operator==(const EdgeEdgeConstraint& other) const
{
    return EdgeEdgeCandidate::operator==(other) && dtype == other.dtype;
}

bool EdgeEdgeConstraint::operator!=(const EdgeEdgeConstraint& other) const
{
    return !(*this == other);
}

bool EdgeEdgeConstraint::operator<(const EdgeEdgeConstraint& other) const
{
    if (EdgeEdgeCandidate::operator==(other)) {
        return dtype < other.dtype;
    }
    return EdgeEdgeCandidate::operator<(other);
}

} // namespace ipc
