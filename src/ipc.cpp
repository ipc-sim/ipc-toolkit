#include "ipc.hpp"

#include <barrier/barrier.hpp>
#include <spatial_hash/hash_grid.hpp>

#include <distance/edge_edge.hpp>
#include <distance/edge_edge_mollifier.hpp>
#include <distance/point_edge.hpp>
#include <distance/point_triangle.hpp>

#include <ccd/edge_vertex_ccd_2D.hpp>
// Etienne Vouga's CCD using a root finder in floating points
#include <CTCD.h>

namespace ipc {

void construct_constraint_set(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    double dhat_squared,
    ccd::Candidates& constraint_set,
    bool ignore_internal_vertices)
{
    double dhat = std::sqrt(dhat_squared);

    ccd::Candidates candidates;
    ccd::HashGrid hash_grid;
    hash_grid.resize(V, V, E, /*inflation_radius=*/dhat);

    // Assumes the edges connect to all boundary vertices
    if (ignore_internal_vertices) {
        for (int e = 0; e < E.rows(); ++e) {
            const int e0 = E(e, 0);
            const int e1 = E(e, 1);
            hash_grid.addVertex(
                V.row(e0), V.row(e0), e0, /*inflation_radius=*/dhat);
            hash_grid.addVertex(
                V.row(e1), V.row(e1), e1, /*inflation_radius=*/dhat);
        }
    } else {
        hash_grid.addVertices(V, V, /*inflation_radius=*/dhat);
    }

    hash_grid.addEdges(V, V, E, /*inflation_radius=*/dhat);
    if (V.cols() == 3) {
        // These are not needed for 2D
        hash_grid.addFaces(V, V, F, /*inflation_radius=*/dhat);
    }

    if (V.cols() == 2) {
        // This is not needed for 3D
        hash_grid.getVertexEdgePairs(
            E, /*group_ids=*/Eigen::VectorXi(), candidates.ev_candidates);
    } else {
        // These are not needed for 2D
        hash_grid.getEdgeEdgePairs(
            E, /*group_ids=*/Eigen::VectorXi(), candidates.ee_candidates);
        hash_grid.getFaceVertexPairs(
            F, /*group_ids=*/Eigen::VectorXi(), candidates.fv_candidates);
    }

    // Cull the candidates by measuring the distance and dropping those that are
    // greater than dhat.

    // TODO: Consider parallelizing this loop
    for (const auto& ev_candidate : candidates.ev_candidates) {
        double distance_sqr = point_edge_distance(
            V.row(ev_candidate.vertex_index),
            V.row(E(ev_candidate.edge_index, 0)),
            V.row(E(ev_candidate.edge_index, 1)));

        if (distance_sqr < dhat_squared) {
            constraint_set.ev_candidates.push_back(ev_candidate);
        }
    }

    for (const auto& ee_candidate : candidates.ee_candidates) {
        double distance_sqr = edge_edge_distance(
            V.row(E(ee_candidate.edge0_index, 0)),
            V.row(E(ee_candidate.edge0_index, 1)),
            V.row(E(ee_candidate.edge1_index, 0)),
            V.row(E(ee_candidate.edge1_index, 1)));

        if (distance_sqr < dhat_squared) {
            constraint_set.ee_candidates.push_back(ee_candidate);
        }
    }

    for (const auto& fv_candidate : candidates.fv_candidates) {
        double distance_sqr = point_triangle_distance(
            V.row(fv_candidate.vertex_index),
            V.row(F(fv_candidate.face_index, 0)),
            V.row(F(fv_candidate.face_index, 1)),
            V.row(F(fv_candidate.face_index, 2)));

        if (distance_sqr < dhat_squared) {
            constraint_set.fv_candidates.push_back(fv_candidate);
        }
    }
}

double compute_barrier_potential(
    const Eigen::MatrixXd& V_rest,
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const ccd::Candidates& constraint_set,
    double dhat_squared)
{
    double potential = 0;

    for (const auto& ev_candidate : constraint_set.ev_candidates) {
        double distance_sqr = point_edge_distance(
            V.row(ev_candidate.vertex_index),
            V.row(E(ev_candidate.edge_index, 0)),
            V.row(E(ev_candidate.edge_index, 1)));
        potential += barrier(distance_sqr, dhat_squared);
    }

    for (const auto& ee_candidate : constraint_set.ee_candidates) {
        const auto& ea0 = V.row(E(ee_candidate.edge0_index, 0));
        const auto& ea1 = V.row(E(ee_candidate.edge0_index, 1));
        const auto& eb0 = V.row(E(ee_candidate.edge1_index, 0));
        const auto& eb1 = V.row(E(ee_candidate.edge1_index, 1));

        double distance_sqr = edge_edge_distance(ea0, ea1, eb0, eb1);
        double eps_x = edge_edge_mollifier_threshold(
            V_rest.row(E(ee_candidate.edge0_index, 0)),
            V_rest.row(E(ee_candidate.edge0_index, 1)),
            V_rest.row(E(ee_candidate.edge1_index, 0)),
            V_rest.row(E(ee_candidate.edge1_index, 1)));
        potential += edge_edge_mollifier(ea0, ea1, eb0, eb1, eps_x)
            * barrier(distance_sqr, dhat_squared);
    }

    for (const auto& fv_candidate : constraint_set.fv_candidates) {
        const auto& p = V.row(fv_candidate.vertex_index);
        const auto& t0 = V.row(F(fv_candidate.face_index, 0));
        const auto& t1 = V.row(F(fv_candidate.face_index, 1));
        const auto& t2 = V.row(F(fv_candidate.face_index, 2));

        double distance_sqr = point_triangle_distance(p, t0, t1, t2);
        potential += barrier(distance_sqr, dhat_squared);
    }

    return potential;
}

Eigen::VectorXd compute_barrier_potential_gradient(
    const Eigen::MatrixXd& V_rest,
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const ccd::Candidates& constraint_set,
    double dhat_squared)
{
    Eigen::VectorXd grad = Eigen::VectorXd::Zero(V.size());
    int dim = V.cols();

    for (const auto& ev_candidate : constraint_set.ev_candidates) {
        // ∇b(d(x)) = b'(d(x)) * ∇d(x)
        const auto& p = V.row(ev_candidate.vertex_index);
        const auto& e0 = V.row(E(ev_candidate.edge_index, 0));
        const auto& e1 = V.row(E(ev_candidate.edge_index, 1));

        Eigen::VectorXd local_grad;
        point_edge_distance_gradient(p, e0, e1, local_grad);

        double distance_sqr = point_edge_distance(p, e0, e1);
        local_grad *= barrier_gradient(distance_sqr, dhat_squared);

        // Map from local to global gradient
        grad.segment(dim * ev_candidate.vertex_index, dim) +=
            local_grad.head(dim);
        grad.segment(dim * E(ev_candidate.edge_index, 0), dim) +=
            local_grad.segment(dim, dim);
        grad.segment(dim * E(ev_candidate.edge_index, 1), dim) +=
            local_grad.tail(dim);
    }

    for (const auto& ee_candidate : constraint_set.ee_candidates) {
        // ∇[m(x) * b(d(x))] = (∇m(x)) * b(d(x)) + m(x) * b'(d(x)) * ∇d(x)
        const auto& ea0 = V.row(E(ee_candidate.edge0_index, 0));
        const auto& ea1 = V.row(E(ee_candidate.edge0_index, 1));
        const auto& eb0 = V.row(E(ee_candidate.edge1_index, 0));
        const auto& eb1 = V.row(E(ee_candidate.edge1_index, 1));

        double distance_sqr = edge_edge_distance(ea0, ea1, eb0, eb1);
        Eigen::VectorXd local_distance_grad;
        edge_edge_distance_gradient(ea0, ea1, eb0, eb1, local_distance_grad);

        double eps_x = edge_edge_mollifier_threshold(
            V_rest.row(E(ee_candidate.edge0_index, 0)),
            V_rest.row(E(ee_candidate.edge0_index, 1)),
            V_rest.row(E(ee_candidate.edge1_index, 0)),
            V_rest.row(E(ee_candidate.edge1_index, 1)));
        double mollifier = edge_edge_mollifier(ea0, ea1, eb0, eb1, eps_x);
        Eigen::VectorXd local_mollifier_grad;
        edge_edge_mollifier_gradient(
            ea0, ea1, eb0, eb1, eps_x, local_mollifier_grad);

        Eigen::VectorXd local_grad =
            local_mollifier_grad * barrier(distance_sqr, dhat_squared)
            + mollifier * barrier_gradient(distance_sqr, dhat_squared)
                * local_distance_grad;

        // Map from local to global gradient
        grad.segment(dim * E(ee_candidate.edge0_index, 0), dim) +=
            local_grad.head(dim);
        grad.segment(dim * E(ee_candidate.edge0_index, 1), dim) +=
            local_grad.segment(dim, dim);
        grad.segment(dim * E(ee_candidate.edge1_index, 0), dim) +=
            local_grad.segment(2 * dim, dim);
        grad.segment(dim * E(ee_candidate.edge1_index, 1), dim) +=
            local_grad.tail(dim);
    }

    for (const auto& fv_candidate : constraint_set.fv_candidates) {
        // ∇b(d(x)) = b'(d(x)) * ∇d(x)
        const auto& p = V.row(fv_candidate.vertex_index);
        const auto& t0 = V.row(F(fv_candidate.face_index, 0));
        const auto& t1 = V.row(F(fv_candidate.face_index, 1));
        const auto& t2 = V.row(F(fv_candidate.face_index, 2));

        Eigen::VectorXd local_grad;
        point_triangle_distance_gradient(p, t0, t1, t2, local_grad);

        double distance_sqr = point_triangle_distance(p, t0, t1, t2);
        local_grad *= barrier_gradient(distance_sqr, dhat_squared);

        // Map from local to global gradient
        grad.segment(dim * fv_candidate.vertex_index, dim) +=
            local_grad.head(dim);
        grad.segment(dim * F(fv_candidate.face_index, 0), dim) +=
            local_grad.segment(dim, dim);
        grad.segment(dim * F(fv_candidate.face_index, 1), dim) +=
            local_grad.segment(2 * dim, dim);
        grad.segment(dim * F(fv_candidate.face_index, 2), dim) +=
            local_grad.tail(dim);
    }

    return grad;
}

void local_hessian_to_global_triplets(
    const Eigen::MatrixXd& local_hessian,
    const std::vector<long>& ids,
    int dim,
    std::vector<Eigen::Triplet<double>>& triplets)
{
    for (int i = 0; i < ids.size(); i++) {
        for (int j = 0; j < ids.size(); j++) {
            for (int k = 0; k < dim; k++) {
                for (int l = 0; l < dim; l++) {
                    triplets.emplace_back(
                        dim * ids[i] + k, dim * ids[j] + l,
                        local_hessian(dim * i + k, dim * j + l));
                }
            }
        }
    }
}

Eigen::SparseMatrix<double> compute_barrier_potential_hessian(
    const Eigen::MatrixXd& V_rest,
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const ccd::Candidates& constraint_set,
    double dhat_squared)
{
    std::vector<Eigen::Triplet<double>> hess_triplets;
    int dim = V.cols();
    hess_triplets.reserve(
        constraint_set.ev_candidates.size() * /*3*3=*/9 * dim * dim
        + constraint_set.ee_candidates.size() * /*4*4=*/16 * dim * dim
        + constraint_set.fv_candidates.size() * /*4*4=*/16 * dim * dim);

    for (const auto& ev_candidate : constraint_set.ev_candidates) {
        // ∇²b(d(x)) = ∇(b'(d(x)) * ∇d(x))
        //           = b"(d(x)) * ∇d(x) * ∇d(x)ᵀ + b'(d(x)) * ∇²d(x)
        const auto& p = V.row(ev_candidate.vertex_index);
        const auto& e0 = V.row(E(ev_candidate.edge_index, 0));
        const auto& e1 = V.row(E(ev_candidate.edge_index, 1));

        double distance_sqr = point_edge_distance(p, e0, e1);
        Eigen::VectorXd local_grad;
        point_edge_distance_gradient(p, e0, e1, local_grad);
        Eigen::MatrixXd local_hess;
        point_edge_distance_hessian(
            p, e0, e1, local_hess, /*project_to_psd=*/true);

        local_hess *= barrier_gradient(distance_sqr, dhat_squared);
        local_hess += barrier_hessian(distance_sqr, dhat_squared) * local_grad
            * local_grad.transpose();

        // Map from local to global gradient
        std::vector<long> ids = { { ev_candidate.vertex_index,
                                    E(ev_candidate.edge_index, 0),
                                    E(ev_candidate.edge_index, 1) } };
        local_hessian_to_global_triplets(local_hess, ids, dim, hess_triplets);
    }

    for (const auto& ee_candidate : constraint_set.ee_candidates) {
        // ∇²[m(x) * b(d(x))] = ∇[∇m(x) * b(d(x)) + m(x) * b'(d(x)) * ∇d(x)]
        //                    = ∇²m(x) * b(d(x)) + b'(d(x)) * ∇d(x) * ∇m(x)ᵀ
        //                      + ∇m(x) * b'(d(x)) * ∇d(x))ᵀ
        //                      + m(x) * b"(d(x)) * ∇d(x) * ∇d(x)ᵀ
        //                      + m(x) * b'(d(x)) * ∇²d(x)
        const auto& ea0 = V.row(E(ee_candidate.edge0_index, 0));
        const auto& ea1 = V.row(E(ee_candidate.edge0_index, 1));
        const auto& eb0 = V.row(E(ee_candidate.edge1_index, 0));
        const auto& eb1 = V.row(E(ee_candidate.edge1_index, 1));

        // Compute distance derivatives
        double distance_sqr = edge_edge_distance(ea0, ea1, eb0, eb1);
        Eigen::VectorXd distance_grad;
        edge_edge_distance_gradient(ea0, ea1, eb0, eb1, distance_grad);
        Eigen::MatrixXd distance_hess;
        edge_edge_distance_hessian(
            ea0, ea1, eb0, eb1, distance_hess, /*project_to_psd=*/true);

        // Compute mollifier derivatives
        double eps_x = edge_edge_mollifier_threshold(
            V_rest.row(E(ee_candidate.edge0_index, 0)),
            V_rest.row(E(ee_candidate.edge0_index, 1)),
            V_rest.row(E(ee_candidate.edge1_index, 0)),
            V_rest.row(E(ee_candidate.edge1_index, 1)));
        double mollifier = edge_edge_mollifier(ea0, ea1, eb0, eb1, eps_x);
        Eigen::VectorXd mollifier_grad;
        edge_edge_mollifier_gradient(ea0, ea1, eb0, eb1, eps_x, mollifier_grad);
        Eigen::MatrixXd mollifier_hess;
        edge_edge_mollifier_hessian(ea0, ea1, eb0, eb1, eps_x, mollifier_hess);

        // Compute_barrier_derivatives
        double b = barrier(distance_sqr, dhat_squared);
        double grad_b = barrier_gradient(distance_sqr, dhat_squared);
        double hess_b = barrier_hessian(distance_sqr, dhat_squared);

        Eigen::MatrixXd local_hess = mollifier_hess * b
            + grad_b
                * (distance_grad * mollifier_grad.transpose()
                   + mollifier_grad * distance_grad.transpose())
            + mollifier
                * (hess_b * distance_grad * distance_grad.transpose()
                   + grad_b * distance_hess);

        // Map from local to global gradient
        std::vector<long> ids = {
            { E(ee_candidate.edge0_index, 0), E(ee_candidate.edge0_index, 1),
              E(ee_candidate.edge1_index, 0), E(ee_candidate.edge1_index, 1) }
        };
        local_hessian_to_global_triplets(local_hess, ids, dim, hess_triplets);
    }

    for (const auto& fv_candidate : constraint_set.fv_candidates) {
        // ∇²b(d(x)) = ∇(b'(d(x)) * ∇d(x))
        //           = b"(d(x)) * ∇d(x) * ∇d(x)ᵀ + b'(d(x)) * ∇²d(x)
        const auto& p = V.row(fv_candidate.vertex_index);
        const auto& t0 = V.row(F(fv_candidate.face_index, 0));
        const auto& t1 = V.row(F(fv_candidate.face_index, 1));
        const auto& t2 = V.row(F(fv_candidate.face_index, 2));

        double distance_sqr = point_triangle_distance(p, t0, t1, t2);
        Eigen::VectorXd local_grad;
        point_triangle_distance_gradient(p, t0, t1, t2, local_grad);
        Eigen::MatrixXd local_hess;
        point_triangle_distance_hessian(
            p, t0, t1, t2, local_hess, /*project_to_psd=*/true);

        local_hess *= barrier_gradient(distance_sqr, dhat_squared);
        local_hess += barrier_hessian(distance_sqr, dhat_squared) * local_grad
            * local_grad.transpose();

        // Map from local to global gradient
        std::vector<long> ids = {
            { fv_candidate.vertex_index, F(fv_candidate.face_index, 0),
              F(fv_candidate.face_index, 1), F(fv_candidate.face_index, 2) }
        };
        local_hessian_to_global_triplets(local_hess, ids, dim, hess_triplets);
    }

    Eigen::SparseMatrix<double> hess(V.size(), V.size());
    hess.setFromTriplets(hess_triplets.begin(), hess_triplets.end());
    return hess;
}

bool is_step_collision_free(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    bool ignore_internal_vertices)
{
    int dim = V0.cols();
    assert(V1.cols() == dim);

    // Broad phase
    ccd::Candidates candidates;
    ccd::HashGrid hash_grid;
    hash_grid.resize(V0, V1, E);

    // Assumes the edges connect to all boundary vertices
    if (ignore_internal_vertices) {
        for (int e = 0; e < E.rows(); ++e) {
            const int e0 = E(e, 0);
            const int e1 = E(e, 1);
            hash_grid.addVertex(V0.row(e0), V1.row(e0), e0);
            hash_grid.addVertex(V0.row(e1), V1.row(e1), e1);
        }
    } else {
        hash_grid.addVertices(V0, V1);
    }
    hash_grid.addEdges(V0, V1, E);
    if (dim == 3) {
        // These are not needed for 2D
        hash_grid.addFaces(V0, V1, F);
    }

    if (dim == 2) {
        // This is not needed for 3D
        hash_grid.getVertexEdgePairs(
            E, /*group_ids=*/Eigen::VectorXi(), candidates.ev_candidates);
    } else {
        // These are not needed for 2D
        hash_grid.getEdgeEdgePairs(
            E, /*group_ids=*/Eigen::VectorXi(), candidates.ee_candidates);
        hash_grid.getFaceVertexPairs(
            F, /*group_ids=*/Eigen::VectorXi(), candidates.fv_candidates);
    }

    // Narrow phase
    double eta = 1e-6;

    for (const auto& ev_candidate : candidates.ev_candidates) {
        double toi;
        double alpha;
        bool is_collision = ccd::compute_edge_vertex_time_of_impact(
            // Displacement of Edge at t=0
            V0.row(E(ev_candidate.edge_index, 0)),
            V0.row(E(ev_candidate.edge_index, 1)),
            // Displacement of Point at t=0
            V0.row(ev_candidate.vertex_index),
            // Displacement of Edge at t=1
            V1.row(E(ev_candidate.edge_index, 0))
                - V0.row(E(ev_candidate.edge_index, 0)),
            V1.row(E(ev_candidate.edge_index, 1))
                - V0.row(E(ev_candidate.edge_index, 1)),
            // Displacement of Point at t=1
            V1.row(ev_candidate.vertex_index)
                - V0.row(ev_candidate.vertex_index),
            toi, alpha);

        if (is_collision) {
            return false;
        }
    }

    for (const auto& ee_candidate : candidates.ee_candidates) {
        double toi;
        bool is_collision = CTCD::edgeEdgeCTCD(
            // Edge 1 at t=0
            V0.row(E(ee_candidate.edge0_index, 0)),
            V0.row(E(ee_candidate.edge0_index, 1)),
            // Edge 2 at t=0
            V0.row(E(ee_candidate.edge1_index, 0)),
            V0.row(E(ee_candidate.edge1_index, 1)),
            // Edge 1 at t=1
            V1.row(E(ee_candidate.edge0_index, 0)),
            V1.row(E(ee_candidate.edge0_index, 1)),
            // Edge 2 at t=1
            V1.row(E(ee_candidate.edge1_index, 0)),
            V1.row(E(ee_candidate.edge1_index, 1)), //
            eta, toi);

        if (is_collision) {
            return false;
        }
    }

    for (const auto& fv_candidate : candidates.fv_candidates) {
        double toi;
        bool is_collision = CTCD::vertexFaceCTCD(
            // Point at t=0
            V0.row(fv_candidate.vertex_index),
            // Triangle at t = 0
            V0.row(F(fv_candidate.face_index, 0)),
            V0.row(F(fv_candidate.face_index, 1)),
            V0.row(F(fv_candidate.face_index, 2)),
            // Point at t=1
            V1.row(fv_candidate.vertex_index),
            // Triangle at t = 1
            V1.row(F(fv_candidate.face_index, 0)),
            V1.row(F(fv_candidate.face_index, 1)),
            V1.row(F(fv_candidate.face_index, 2)), //
            eta, toi);

        if (is_collision) {
            return false;
        }
    }

    return true;
}

double compute_collision_free_stepsize(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    bool ignore_internal_vertices)
{
    int dim = V0.cols();
    assert(V1.cols() == dim);

    // Broad phase
    ccd::Candidates candidates;
    ccd::HashGrid hash_grid;
    hash_grid.resize(V0, V1, E);

    // Assumes the edges connect to all boundary vertices
    if (ignore_internal_vertices) {
        for (int e = 0; e < E.rows(); ++e) {
            const int e0 = E(e, 0);
            const int e1 = E(e, 1);
            hash_grid.addVertex(V0.row(e0), V1.row(e0), e0);
            hash_grid.addVertex(V0.row(e1), V1.row(e1), e1);
        }
    } else {
        hash_grid.addVertices(V0, V1);
    }
    hash_grid.addEdges(V0, V1, E);
    if (dim == 3) {
        // These are not needed for 2D
        hash_grid.addFaces(V0, V1, F);
    }

    if (dim == 2) {
        // This is not needed for 3D
        hash_grid.getVertexEdgePairs(
            E, /*group_ids=*/Eigen::VectorXi(), candidates.ev_candidates);
    } else {
        // These are not needed for 2D
        hash_grid.getEdgeEdgePairs(
            E, /*group_ids=*/Eigen::VectorXi(), candidates.ee_candidates);
        hash_grid.getFaceVertexPairs(
            F, /*group_ids=*/Eigen::VectorXi(), candidates.fv_candidates);
    }

    // Narrow phase
    const double eta = 1e-6;
    double earliest_toi = std::numeric_limits<double>::infinity();

    for (const auto& ev_candidate : candidates.ev_candidates) {
        double toi;
        double alpha;
        bool is_collision = ccd::compute_edge_vertex_time_of_impact(
            // Edge at t=0
            V0.row(E(ev_candidate.edge_index, 0)),
            V0.row(E(ev_candidate.edge_index, 1)),
            // Point at t=0
            V0.row(ev_candidate.vertex_index),
            // Displacement of Edge at t=1
            V1.row(E(ev_candidate.edge_index, 0))
                - V0.row(E(ev_candidate.edge_index, 0)),
            V1.row(E(ev_candidate.edge_index, 1))
                - V0.row(E(ev_candidate.edge_index, 1)), //
            // Displacement of Point at t=1
            V1.row(ev_candidate.vertex_index)
                - V0.row(ev_candidate.vertex_index),
            toi, alpha);

        if (is_collision && toi < earliest_toi) {
            earliest_toi = toi;
        }
    }

    for (const auto& ee_candidate : candidates.ee_candidates) {
        double toi;
        bool is_collision = CTCD::edgeEdgeCTCD(
            // Edge 1 at t=0
            V0.row(E(ee_candidate.edge0_index, 0)),
            V0.row(E(ee_candidate.edge0_index, 1)),
            // Edge 2 at t=0
            V0.row(E(ee_candidate.edge1_index, 0)),
            V0.row(E(ee_candidate.edge1_index, 1)),
            // Edge 1 at t=1
            V1.row(E(ee_candidate.edge0_index, 0)),
            V1.row(E(ee_candidate.edge0_index, 1)),
            // Edge 2 at t=1
            V1.row(E(ee_candidate.edge1_index, 0)),
            V1.row(E(ee_candidate.edge1_index, 1)), //
            eta, toi);

        if (is_collision && toi < earliest_toi) {
            earliest_toi = toi;
        }
    }

    for (const auto& fv_candidate : candidates.fv_candidates) {
        double toi;
        bool is_collision = CTCD::vertexFaceCTCD(
            // Point at t=0
            V0.row(fv_candidate.vertex_index),
            // Triangle at t = 0
            V0.row(F(fv_candidate.face_index, 0)),
            V0.row(F(fv_candidate.face_index, 1)),
            V0.row(F(fv_candidate.face_index, 2)),
            // Point at t=1
            V1.row(fv_candidate.vertex_index),
            // Triangle at t = 1
            V1.row(F(fv_candidate.face_index, 0)),
            V1.row(F(fv_candidate.face_index, 1)),
            V1.row(F(fv_candidate.face_index, 2)), //
            eta, toi);

        if (is_collision && toi < earliest_toi) {
            earliest_toi = toi;
        }
    }

    return earliest_toi;
}

// NOTE: Actually distance squared
double compute_minimum_distance(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const ccd::Candidates& constraint_set)
{
    double min_distance = std::numeric_limits<double>::infinity();

    for (const auto& ev_candidate : constraint_set.ev_candidates) {
        double distance_sqr = point_edge_distance(
            V.row(ev_candidate.vertex_index),
            V.row(E(ev_candidate.edge_index, 0)),
            V.row(E(ev_candidate.edge_index, 1)));

        if (distance_sqr < min_distance) {
            min_distance = distance_sqr;
        }
    }

    for (const auto& ee_candidate : constraint_set.ee_candidates) {
        double distance_sqr = edge_edge_distance(
            V.row(E(ee_candidate.edge0_index, 0)),
            V.row(E(ee_candidate.edge0_index, 1)),
            V.row(E(ee_candidate.edge1_index, 0)),
            V.row(E(ee_candidate.edge1_index, 1)));

        if (distance_sqr < min_distance) {
            min_distance = distance_sqr;
        }
    }

    for (const auto& fv_candidate : constraint_set.fv_candidates) {
        double distance_sqr = point_triangle_distance(
            V.row(fv_candidate.vertex_index),
            V.row(F(fv_candidate.face_index, 0)),
            V.row(F(fv_candidate.face_index, 1)),
            V.row(F(fv_candidate.face_index, 2)));

        if (distance_sqr < min_distance) {
            min_distance = distance_sqr;
        }
    }

    return min_distance;
}

} // namespace ipc
