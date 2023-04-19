#include "edge_edge.hpp"

#include <ipc/barrier/barrier.hpp>
#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/edge_edge_mollifier.hpp>

namespace ipc {

EdgeEdgeConstraint::EdgeEdgeConstraint(
    long edge0_id, long edge1_id, double eps_x)
    : EdgeEdgeCandidate(edge0_id, edge1_id)
    , eps_x(eps_x)
{
}

EdgeEdgeConstraint::EdgeEdgeConstraint(
    const EdgeEdgeCandidate& candidate, double eps_x)
    : EdgeEdgeCandidate(candidate)
    , eps_x(eps_x)
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
    const double adjusted_dhat = 2 * minimum_distance * dhat + dhat * dhat;
    const double min_dist_squared = minimum_distance * minimum_distance;

    // ∇[m(x) * b(d(x))] = (∇m(x)) * b(d(x)) + m(x) * b'(d(x)) * ∇d(x)
    const auto& [ea0, ea1, eb0, eb1] = this->vertices(vertices, edges, faces);

    // The distance type is unknown because of mollified PP and PE
    // constraints where also added as EE constraints.
    const EdgeEdgeDistanceType dtype =
        edge_edge_distance_type(ea0, ea1, eb0, eb1);
    const double distance = edge_edge_distance(ea0, ea1, eb0, eb1, dtype);
    const Vector12d distance_grad =
        edge_edge_distance_gradient(ea0, ea1, eb0, eb1, dtype);

    // m(x)
    const double mollifier = edge_edge_mollifier(ea0, ea1, eb0, eb1, eps_x);
    // ∇m(x)
    const Vector12d mollifier_grad =
        edge_edge_mollifier_gradient(ea0, ea1, eb0, eb1, eps_x);

    // b(d(x))
    const double b = barrier(distance - min_dist_squared, adjusted_dhat);
    // b'(d(x))
    const double grad_b =
        barrier_gradient(distance - min_dist_squared, adjusted_dhat);

    return weight * (mollifier_grad * b + mollifier * grad_b * distance_grad);
}

MatrixMax12d EdgeEdgeConstraint::compute_potential_hessian(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const double dhat,
    const bool project_hessian_to_psd) const
{
    const double adjusted_dhat = 2 * minimum_distance * dhat + dhat * dhat;
    const double min_dist_squared = minimum_distance * minimum_distance;

    // ∇²[m(x) * b(d(x))] = ∇[∇m(x) * b(d(x)) + m(x) * b'(d(x)) * ∇d(x)]
    //                    = ∇²m(x) * b(d(x)) + b'(d(x)) * ∇d(x) * ∇m(x)ᵀ
    //                      + ∇m(x) * b'(d(x)) * ∇d(x))ᵀ
    //                      + m(x) * b"(d(x)) * ∇d(x) * ∇d(x)ᵀ
    //                      + m(x) * b'(d(x)) * ∇²d(x)
    const auto& [ea0, ea1, eb0, eb1] = this->vertices(vertices, edges, faces);

    // Compute distance derivatives
    // The distance type is unknown because of mollified PP and PE
    // constraints where also added as EE constraints.
    const EdgeEdgeDistanceType dtype =
        edge_edge_distance_type(ea0, ea1, eb0, eb1);
    const double distance = edge_edge_distance(ea0, ea1, eb0, eb1, dtype);
    const Vector12d distance_grad =
        edge_edge_distance_gradient(ea0, ea1, eb0, eb1, dtype);
    const Matrix12d distance_hess =
        edge_edge_distance_hessian(ea0, ea1, eb0, eb1, dtype);

    // Compute mollifier derivatives
    const double mollifier = edge_edge_mollifier(ea0, ea1, eb0, eb1, eps_x);
    const VectorMax12d mollifier_grad =
        edge_edge_mollifier_gradient(ea0, ea1, eb0, eb1, eps_x);
    const MatrixMax12d mollifier_hess =
        edge_edge_mollifier_hessian(ea0, ea1, eb0, eb1, eps_x);

    // Compute barrier derivatives
    const double b = barrier(distance - min_dist_squared, adjusted_dhat);
    const double grad_b =
        barrier_gradient(distance - min_dist_squared, adjusted_dhat);
    const double hess_b =
        barrier_hessian(distance - min_dist_squared, adjusted_dhat);

    MatrixMax12d hess = mollifier_hess * b
        + grad_b
            * (distance_grad * mollifier_grad.transpose()
               + mollifier_grad * distance_grad.transpose())
        + mollifier
            * (hess_b * distance_grad * distance_grad.transpose()
               + grad_b * distance_hess);

    if (project_hessian_to_psd) {
        hess = project_to_psd(hess);
    }

    return weight * hess;
}

} // namespace ipc
