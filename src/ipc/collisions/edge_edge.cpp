#include "edge_edge.hpp"

#include <ipc/barrier/barrier.hpp>
#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/edge_edge_mollifier.hpp>

namespace ipc {

EdgeEdgeConstraint::EdgeEdgeConstraint(
    long edge0_index, long edge1_index, double eps_x)
    : EdgeEdgeCandidate(edge0_index, edge1_index)
    , eps_x(eps_x)
{
}

EdgeEdgeConstraint::EdgeEdgeConstraint(
    const EdgeEdgeCandidate& candidate, double eps_x)
    : EdgeEdgeCandidate(candidate)
    , eps_x(eps_x)
{
}

double EdgeEdgeConstraint::compute_distance(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F) const
{
    // The distance type is unknown because of mollified PP and PE
    // constraints where also added as EE constraints.
    return edge_edge_distance(
        V.row(E(edge0_index, 0)), V.row(E(edge0_index, 1)),
        V.row(E(edge1_index, 0)), V.row(E(edge1_index, 1)));
}

VectorMax12d EdgeEdgeConstraint::compute_distance_gradient(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F) const
{
    VectorMax12d distance_grad;
    edge_edge_distance_gradient(
        V.row(E(edge0_index, 0)), V.row(E(edge0_index, 1)),
        V.row(E(edge1_index, 0)), V.row(E(edge1_index, 1)), distance_grad);
    return distance_grad;
}

MatrixMax12d EdgeEdgeConstraint::compute_distance_hessian(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F) const
{
    MatrixMax12d distance_hess;
    edge_edge_distance_hessian(
        V.row(E(edge0_index, 0)), V.row(E(edge0_index, 1)),
        V.row(E(edge1_index, 0)), V.row(E(edge1_index, 1)), distance_hess);
    return distance_hess;
}

double EdgeEdgeConstraint::compute_potential(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const double dhat) const
{
    return edge_edge_mollifier(
               V.row(E(edge0_index, 0)), V.row(E(edge0_index, 1)),
               V.row(E(edge1_index, 0)), V.row(E(edge1_index, 1)), eps_x)
        * CollisionConstraint::compute_potential(V, E, F, dhat);
}

VectorMax12d EdgeEdgeConstraint::compute_potential_gradient(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const double dhat) const
{
    const double dhat_squared = dhat * dhat;

    // ∇[m(x) * b(d(x))] = (∇m(x)) * b(d(x)) + m(x) * b'(d(x)) * ∇d(x)
    const auto& ea0 = V.row(E(edge0_index, 0));
    const auto& ea1 = V.row(E(edge0_index, 1));
    const auto& eb0 = V.row(E(edge1_index, 0));
    const auto& eb1 = V.row(E(edge1_index, 1));

    // The distance type is unknown because of mollified PP and PE
    // constraints where also added as EE constraints.
    const EdgeEdgeDistanceType dtype =
        edge_edge_distance_type(ea0, ea1, eb0, eb1);
    const double distance = edge_edge_distance(ea0, ea1, eb0, eb1, dtype);
    VectorMax12d distance_grad;
    edge_edge_distance_gradient(ea0, ea1, eb0, eb1, dtype, distance_grad);

    // m(x)
    const double mollifier = edge_edge_mollifier(ea0, ea1, eb0, eb1, eps_x);
    // ∇m(x)
    VectorMax12d mollifier_grad;
    edge_edge_mollifier_gradient(ea0, ea1, eb0, eb1, eps_x, mollifier_grad);

    // b(d(x))
    const double b = barrier(
        distance - minimum_distance * minimum_distance,
        2 * minimum_distance * dhat + dhat_squared);
    // b'(d(x))
    const double grad_b = barrier_gradient(
        distance - minimum_distance * minimum_distance,
        2 * minimum_distance * dhat + dhat_squared);

    return weight * (mollifier_grad * b + mollifier * grad_b * distance_grad);
}

MatrixMax12d EdgeEdgeConstraint::compute_potential_hessian(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const double dhat,
    const bool project_hessian_to_psd) const
{
    const double dhat_squared = dhat * dhat;
    const double min_dist_squrared = minimum_distance * minimum_distance;

    // ∇²[m(x) * b(d(x))] = ∇[∇m(x) * b(d(x)) + m(x) * b'(d(x)) * ∇d(x)]
    //                    = ∇²m(x) * b(d(x)) + b'(d(x)) * ∇d(x) * ∇m(x)ᵀ
    //                      + ∇m(x) * b'(d(x)) * ∇d(x))ᵀ
    //                      + m(x) * b"(d(x)) * ∇d(x) * ∇d(x)ᵀ
    //                      + m(x) * b'(d(x)) * ∇²d(x)
    const auto& ea0 = V.row(E(edge0_index, 0));
    const auto& ea1 = V.row(E(edge0_index, 1));
    const auto& eb0 = V.row(E(edge1_index, 0));
    const auto& eb1 = V.row(E(edge1_index, 1));

    // Compute distance derivatives
    // The distance type is unknown because of mollified PP and PE
    // constraints where also added as EE constraints.
    const EdgeEdgeDistanceType dtype =
        edge_edge_distance_type(ea0, ea1, eb0, eb1);
    const double distance = edge_edge_distance(ea0, ea1, eb0, eb1, dtype);
    VectorMax12d distance_grad;
    edge_edge_distance_gradient(ea0, ea1, eb0, eb1, dtype, distance_grad);
    MatrixMax12d distance_hess;
    edge_edge_distance_hessian(ea0, ea1, eb0, eb1, dtype, distance_hess);

    // Compute mollifier derivatives
    const double mollifier = edge_edge_mollifier(ea0, ea1, eb0, eb1, eps_x);
    VectorMax12d mollifier_grad;
    edge_edge_mollifier_gradient(ea0, ea1, eb0, eb1, eps_x, mollifier_grad);
    MatrixMax12d mollifier_hess;
    edge_edge_mollifier_hessian(ea0, ea1, eb0, eb1, eps_x, mollifier_hess);

    // Compute barrier derivatives
    const double b = barrier(
        distance - min_dist_squrared,
        2 * minimum_distance * dhat + dhat_squared);
    const double grad_b = barrier_gradient(
        distance - min_dist_squrared,
        2 * minimum_distance * dhat + dhat_squared);
    const double hess_b = barrier_hessian(
        distance - min_dist_squrared,
        2 * minimum_distance * dhat + dhat_squared);

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
