#include <ipc/ipc.hpp>

#include <unordered_map>

#include <ipc/barrier/barrier.hpp>
#include <ipc/ccd/ccd.hpp>
#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/edge_edge_mollifier.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/distance/point_triangle.hpp>
#include <ipc/spatial_hash/hash_grid.hpp>
#include <ipc/utils/local_hessian_to_global_triplets.hpp>

namespace ipc {

/// Find the index of the undirected edge (e0, e1)
long find_edge(const Eigen::MatrixXi& E, long e0, long e1)
{
    for (long i = 0; i < E.rows(); i++) {
        if ((E(i, 0) == e0 && E(i, 1) == e1)
            || (E(i, 1) == e0 && E(i, 0) == e1)) {
            return i;
        }
    }
    throw "Unable to find edge!";
}

template <typename Hash>
void add_vertex_vertex_constraint(
    std::vector<VertexVertexConstraint>& vv_constraints,
    std::unordered_map<VertexVertexConstraint, long, Hash>& vv_to_index,
    const long v0i,
    const long v1i)
{
    VertexVertexConstraint vv_constraint(v0i, v1i);
    auto found_item = vv_to_index.find(vv_constraint);
    if (found_item != vv_to_index.end()) {
        // Constraint already exists, so increase multiplicity
        vv_constraints[found_item->second].multiplicity++;
    } else {
        // New constraint, so add it to the end of vv_constraints
        vv_to_index.emplace(vv_constraint, vv_constraints.size());
        vv_constraints.push_back(vv_constraint);
    }
}

template <typename Hash>
void add_edge_vertex_constraint(
    std::vector<EdgeVertexConstraint>& ev_constraints,
    std::unordered_map<EdgeVertexConstraint, long, Hash>& ev_to_index,
    const long ei,
    const long vi)
{
    EdgeVertexConstraint ev_constraint(ei, vi);
    auto found_item = ev_to_index.find(ev_constraint);
    if (found_item != ev_to_index.end()) {
        // Constraint already exists, so increase multiplicity
        ev_constraints[found_item->second].multiplicity++;
    } else {
        // New constraint, so add it to the end of vv_constraints
        ev_to_index.emplace(ev_constraint, ev_constraints.size());
        ev_constraints.push_back(ev_constraint);
    }
}

void construct_constraint_set(
    const Eigen::MatrixXd& V_rest,
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    double dhat,
    Constraints& constraint_set,
    bool ignore_internal_vertices,
    const Eigen::VectorXi& vertex_group_ids)
{
    double dhat_squared = dhat * dhat;

    Candidates candidates;
    HashGrid hash_grid;
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
            E, vertex_group_ids, candidates.ev_candidates);
    } else {
        // These are not needed for 2D
        hash_grid.getEdgeEdgePairs(
            E, vertex_group_ids, candidates.ee_candidates);
        hash_grid.getFaceVertexPairs(
            F, vertex_group_ids, candidates.fv_candidates);
    }

    // Cull the candidates by measuring the distance and dropping those that are
    // greater than dhat.

    // Store the indices to VV and EV pairs to avoid duplicates.
    auto vv_hash = [&V](const VertexVertexConstraint& vv_constraint) -> size_t {
        // The vertex-vertex pair should be order independent
        long min_vi =
            std::min(vv_constraint.vertex0_index, vv_constraint.vertex1_index);
        long max_vi =
            std::max(vv_constraint.vertex0_index, vv_constraint.vertex1_index);
        return size_t(max_vi * V.rows() + min_vi);
    };
    std::unordered_map<VertexVertexConstraint, long, decltype(vv_hash)>
        vv_to_index(/*min_buckets=*/candidates.size(), vv_hash);
    auto ev_hash = [&V](const EdgeVertexConstraint& ev_constraint) -> size_t {
        // There are max E.rows() * V.rows() constraints
        return size_t(
            ev_constraint.edge_index * V.rows() + ev_constraint.vertex_index);
    };
    std::unordered_map<EdgeVertexConstraint, long, decltype(ev_hash)>
        ev_to_index(/*min_buckets=*/candidates.size(), ev_hash);

    for (const auto& ev_candidate : candidates.ev_candidates) {
        long vi = ev_candidate.vertex_index;
        long e0i = E(ev_candidate.edge_index, 0);
        long e1i = E(ev_candidate.edge_index, 1);

        PointEdgeDistanceType dtype =
            point_edge_distance_type(V.row(vi), V.row(e0i), V.row(e1i));

        double distance_sqr =
            point_edge_distance(V.row(vi), V.row(e0i), V.row(e1i), dtype);

        if (distance_sqr < dhat_squared) {
            switch (dtype) {
            case PointEdgeDistanceType::P_E0:
                add_vertex_vertex_constraint(
                    constraint_set.vv_constraints, vv_to_index, vi, e0i);
                break;

            case PointEdgeDistanceType::P_E1:
                add_vertex_vertex_constraint(
                    constraint_set.vv_constraints, vv_to_index, vi, e1i);
                break;

            case PointEdgeDistanceType::P_E:
                // ev_candidates is a set, so no duplicate EV constraints
                ev_to_index.emplace(
                    EdgeVertexCandidate(ev_candidate.edge_index, vi),
                    constraint_set.ev_constraints.size());
                constraint_set.ev_constraints.emplace_back(ev_candidate);
                break;
            }
        }
    }

    for (const auto& ee_candidate : candidates.ee_candidates) {
        long eai = ee_candidate.edge0_index, ebi = ee_candidate.edge1_index;
        long ea0i = E(eai, 0), ea1i = E(eai, 1);
        long eb0i = E(ebi, 0), eb1i = E(ebi, 1);

        EdgeEdgeDistanceType dtype = edge_edge_distance_type(
            V.row(ea0i), V.row(ea1i), V.row(eb0i), V.row(eb1i));

        double distance_sqr = edge_edge_distance(
            V.row(ea0i), V.row(ea1i), V.row(eb0i), V.row(eb1i), dtype);

        if (distance_sqr < dhat_squared) {
            double eps_x = edge_edge_mollifier_threshold(
                V_rest.row(ea0i), V_rest.row(ea1i), //
                V_rest.row(eb0i), V_rest.row(eb1i));
            double ee_cross_norm_sqr = edge_edge_cross_squarednorm(
                V.row(ea0i), V.row(ea1i), V.row(eb0i), V.row(eb1i));
            if (ee_cross_norm_sqr < eps_x) {
                // NOTE: This may not actually the distance type, but all EE
                // pairs requiring mollification must be mollified later.
                dtype = EdgeEdgeDistanceType::EA_EB;
            }

            switch (dtype) {
            case EdgeEdgeDistanceType::EA0_EB0:
                add_vertex_vertex_constraint(
                    constraint_set.vv_constraints, vv_to_index, ea0i, eb0i);
                break;

            case EdgeEdgeDistanceType::EA0_EB1:
                add_vertex_vertex_constraint(
                    constraint_set.vv_constraints, vv_to_index, ea0i, eb1i);
                break;

            case EdgeEdgeDistanceType::EA1_EB0:
                add_vertex_vertex_constraint(
                    constraint_set.vv_constraints, vv_to_index, ea1i, eb0i);
                break;

            case EdgeEdgeDistanceType::EA1_EB1:
                add_vertex_vertex_constraint(
                    constraint_set.vv_constraints, vv_to_index, ea1i, eb1i);
                break;

            case EdgeEdgeDistanceType::EA_EB0:
                add_edge_vertex_constraint(
                    constraint_set.ev_constraints, ev_to_index, eai, eb0i);
                break;

            case EdgeEdgeDistanceType::EA_EB1:
                add_edge_vertex_constraint(
                    constraint_set.ev_constraints, ev_to_index, eai, eb1i);
                break;

            case EdgeEdgeDistanceType::EA0_EB:
                add_edge_vertex_constraint(
                    constraint_set.ev_constraints, ev_to_index, ebi, ea0i);
                break;

            case EdgeEdgeDistanceType::EA1_EB:
                add_edge_vertex_constraint(
                    constraint_set.ev_constraints, ev_to_index, ebi, ea1i);
                break;

            case EdgeEdgeDistanceType::EA_EB:
                constraint_set.ee_constraints.emplace_back(ee_candidate, eps_x);
                break;
            }
        }
    }

    for (const auto& fv_candidate : candidates.fv_candidates) {
        long f0i = F(fv_candidate.face_index, 0);
        long f1i = F(fv_candidate.face_index, 1);
        long f2i = F(fv_candidate.face_index, 2);

        // Compute distance type
        PointTriangleDistanceType dtype = point_triangle_distance_type(
            V.row(fv_candidate.vertex_index), //
            V.row(f0i), V.row(f1i), V.row(f2i));

        double distance_sqr = point_triangle_distance(
            V.row(fv_candidate.vertex_index), //
            V.row(f0i), V.row(f1i), V.row(f2i), dtype);

        if (distance_sqr < dhat_squared) {
            switch (dtype) {
            case PointTriangleDistanceType::P_T0:
                constraint_set.vv_constraints.emplace_back(
                    fv_candidate.vertex_index, f0i);
                break;

            case PointTriangleDistanceType::P_T1:
                constraint_set.vv_constraints.emplace_back(
                    fv_candidate.vertex_index, f1i);
                break;

            case PointTriangleDistanceType::P_T2:
                constraint_set.vv_constraints.emplace_back(
                    fv_candidate.vertex_index, f2i);
                break;

            case PointTriangleDistanceType::P_E0:
                constraint_set.ev_constraints.emplace_back(
                    find_edge(E, f0i, f1i), fv_candidate.vertex_index);
                break;

            case PointTriangleDistanceType::P_E1:
                constraint_set.ev_constraints.emplace_back(
                    find_edge(E, f1i, f2i), fv_candidate.vertex_index);
                break;

            case PointTriangleDistanceType::P_E2:
                constraint_set.ev_constraints.emplace_back(
                    find_edge(E, f2i, f0i), fv_candidate.vertex_index);
                break;

            case PointTriangleDistanceType::P_T:
                constraint_set.fv_constraints.emplace_back(fv_candidate);
                break;
            }
        }
    }
}

double compute_barrier_potential(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const Constraints& constraint_set,
    double dhat)
{
    double dhat_squared = dhat * dhat;
    double potential = 0;

    for (const auto& vv_constraint : constraint_set.vv_constraints) {
        double distance_sqr = point_point_distance(
            V.row(vv_constraint.vertex0_index),
            V.row(vv_constraint.vertex1_index));
        potential +=
            vv_constraint.multiplicity * barrier(distance_sqr, dhat_squared);
    }

    for (const auto& ev_constraint : constraint_set.ev_constraints) {
        // The distance type is known because of construct_constraint_set()
        double distance_sqr = point_edge_distance(
            V.row(ev_constraint.vertex_index),
            V.row(E(ev_constraint.edge_index, 0)),
            V.row(E(ev_constraint.edge_index, 1)), PointEdgeDistanceType::P_E);
        potential +=
            ev_constraint.multiplicity * barrier(distance_sqr, dhat_squared);
    }

    for (const auto& ee_constraint : constraint_set.ee_constraints) {
        const auto& ea0 = V.row(E(ee_constraint.edge0_index, 0));
        const auto& ea1 = V.row(E(ee_constraint.edge0_index, 1));
        const auto& eb0 = V.row(E(ee_constraint.edge1_index, 0));
        const auto& eb1 = V.row(E(ee_constraint.edge1_index, 1));

        // The distance type is unknown because of mollified PP and PE
        // constraints where also added as EE constraints.
        double distance_sqr = edge_edge_distance(ea0, ea1, eb0, eb1);
        potential +=
            edge_edge_mollifier(ea0, ea1, eb0, eb1, ee_constraint.eps_x)
            * barrier(distance_sqr, dhat_squared);
    }

    for (const auto& fv_constraint : constraint_set.fv_constraints) {
        const auto& p = V.row(fv_constraint.vertex_index);
        const auto& t0 = V.row(F(fv_constraint.face_index, 0));
        const auto& t1 = V.row(F(fv_constraint.face_index, 1));
        const auto& t2 = V.row(F(fv_constraint.face_index, 2));

        // The distance type is known because of construct_constraint_set()
        double distance_sqr = point_triangle_distance(
            p, t0, t1, t2, PointTriangleDistanceType::P_T);
        potential += barrier(distance_sqr, dhat_squared);
    }

    return potential;
}

Eigen::VectorXd compute_barrier_potential_gradient(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const Constraints& constraint_set,
    double dhat)
{
    double dhat_squared = dhat * dhat;

    Eigen::VectorXd grad = Eigen::VectorXd::Zero(V.size());
    int dim = V.cols();

    for (const auto& vv_constraint : constraint_set.vv_constraints) {
        // ∇[m * b(d(x))] = m * b'(d(x)) * ∇d(x)
        const auto& p0 = V.row(vv_constraint.vertex0_index);
        const auto& p1 = V.row(vv_constraint.vertex1_index);

        Eigen::VectorXd local_grad;
        point_point_distance_gradient(p0, p1, local_grad);

        double distance_sqr = point_point_distance(p0, p1);
        local_grad *= barrier_gradient(distance_sqr, dhat_squared);

        local_grad *= vv_constraint.multiplicity;

        // Map from local to global gradient
        grad.segment(dim * vv_constraint.vertex0_index, dim) +=
            local_grad.head(dim);
        grad.segment(dim * vv_constraint.vertex1_index, dim) +=
            local_grad.tail(dim);
    }

    for (const auto& ev_constraint : constraint_set.ev_constraints) {
        // ∇[m * b(d(x))] = m * b'(d(x)) * ∇d(x)
        const auto& p = V.row(ev_constraint.vertex_index);
        const auto& e0 = V.row(E(ev_constraint.edge_index, 0));
        const auto& e1 = V.row(E(ev_constraint.edge_index, 1));

        Eigen::VectorXd local_grad;
        point_edge_distance_gradient(
            p, e0, e1, PointEdgeDistanceType::P_E, local_grad);

        double distance_sqr =
            point_edge_distance(p, e0, e1, PointEdgeDistanceType::P_E);
        local_grad *= barrier_gradient(distance_sqr, dhat_squared);

        local_grad *= ev_constraint.multiplicity;

        // Map from local to global gradient
        grad.segment(dim * ev_constraint.vertex_index, dim) +=
            local_grad.head(dim);
        grad.segment(dim * E(ev_constraint.edge_index, 0), dim) +=
            local_grad.segment(dim, dim);
        grad.segment(dim * E(ev_constraint.edge_index, 1), dim) +=
            local_grad.tail(dim);
    }

    for (const auto& ee_constraint : constraint_set.ee_constraints) {
        // ∇[m(x) * b(d(x))] = (∇m(x)) * b(d(x)) + m(x) * b'(d(x)) * ∇d(x)
        const auto& ea0 = V.row(E(ee_constraint.edge0_index, 0));
        const auto& ea1 = V.row(E(ee_constraint.edge0_index, 1));
        const auto& eb0 = V.row(E(ee_constraint.edge1_index, 0));
        const auto& eb1 = V.row(E(ee_constraint.edge1_index, 1));

        // The distance type is unknown because of mollified PP and PE
        // constraints where also added as EE constraints.
        EdgeEdgeDistanceType dtype =
            edge_edge_distance_type(ea0, ea1, eb0, eb1);
        double distance_sqr = edge_edge_distance(ea0, ea1, eb0, eb1, dtype);
        Eigen::VectorXd local_distance_grad;
        edge_edge_distance_gradient(
            ea0, ea1, eb0, eb1, dtype, local_distance_grad);

        double mollifier =
            edge_edge_mollifier(ea0, ea1, eb0, eb1, ee_constraint.eps_x);
        Eigen::VectorXd local_mollifier_grad;
        edge_edge_mollifier_gradient(
            ea0, ea1, eb0, eb1, ee_constraint.eps_x, local_mollifier_grad);

        Eigen::VectorXd local_grad =
            local_mollifier_grad * barrier(distance_sqr, dhat_squared)
            + mollifier * barrier_gradient(distance_sqr, dhat_squared)
                * local_distance_grad;

        // Map from local to global gradient
        grad.segment(dim * E(ee_constraint.edge0_index, 0), dim) +=
            local_grad.head(dim);
        grad.segment(dim * E(ee_constraint.edge0_index, 1), dim) +=
            local_grad.segment(dim, dim);
        grad.segment(dim * E(ee_constraint.edge1_index, 0), dim) +=
            local_grad.segment(2 * dim, dim);
        grad.segment(dim * E(ee_constraint.edge1_index, 1), dim) +=
            local_grad.tail(dim);
    }

    for (const auto& fv_constraint : constraint_set.fv_constraints) {
        // ∇b(d(x)) = b'(d(x)) * ∇d(x)
        const auto& p = V.row(fv_constraint.vertex_index);
        const auto& t0 = V.row(F(fv_constraint.face_index, 0));
        const auto& t1 = V.row(F(fv_constraint.face_index, 1));
        const auto& t2 = V.row(F(fv_constraint.face_index, 2));

        Eigen::VectorXd local_grad;
        point_triangle_distance_gradient(
            p, t0, t1, t2, PointTriangleDistanceType::P_T, local_grad);

        double distance_sqr = point_triangle_distance(
            p, t0, t1, t2, PointTriangleDistanceType::P_T);
        local_grad *= barrier_gradient(distance_sqr, dhat_squared);

        // Map from local to global gradient
        grad.segment(dim * fv_constraint.vertex_index, dim) +=
            local_grad.head(dim);
        grad.segment(dim * F(fv_constraint.face_index, 0), dim) +=
            local_grad.segment(dim, dim);
        grad.segment(dim * F(fv_constraint.face_index, 1), dim) +=
            local_grad.segment(2 * dim, dim);
        grad.segment(dim * F(fv_constraint.face_index, 2), dim) +=
            local_grad.tail(dim);
    }

    return grad;
}

Eigen::SparseMatrix<double> compute_barrier_potential_hessian(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const Constraints& constraint_set,
    double dhat,
    bool project_to_psd)
{
    double dhat_squared = dhat * dhat;

    std::vector<Eigen::Triplet<double>> hess_triplets;
    int dim = V.cols();
    hess_triplets.reserve(
        constraint_set.vv_constraints.size() * /*2*2=*/4 * dim * dim
        + constraint_set.ev_constraints.size() * /*3*3=*/9 * dim * dim
        + constraint_set.ee_constraints.size() * /*4*4=*/16 * dim * dim
        + constraint_set.fv_constraints.size() * /*4*4=*/16 * dim * dim);

    for (const auto& vv_constraint : constraint_set.vv_constraints) {
        // ∇²[m * b(d(x))] = m * ∇(b'(d(x)) * ∇d(x))
        //                 = m * [b"(d(x)) * ∇d(x) * ∇d(x)ᵀ + b'(d(x)) * ∇²d(x)]
        const auto& p0 = V.row(vv_constraint.vertex0_index);
        const auto& p1 = V.row(vv_constraint.vertex1_index);

        double distance_sqr = point_point_distance(p0, p1);
        Eigen::VectorXd local_grad;
        point_point_distance_gradient(p0, p1, local_grad);
        Eigen::MatrixXd local_hess;
        point_point_distance_hessian(p0, p1, local_hess);

        local_hess *= barrier_gradient(distance_sqr, dhat_squared);
        local_hess += barrier_hessian(distance_sqr, dhat_squared) * local_grad
            * local_grad.transpose();

        local_hess *= vv_constraint.multiplicity;

        if (project_to_psd) {
            local_hess = Eigen::project_to_psd(local_hess);
        }

        // Map from local to global gradient
        std::vector<long> ids = { { vv_constraint.vertex0_index,
                                    vv_constraint.vertex1_index } };
        local_hessian_to_global_triplets(local_hess, ids, dim, hess_triplets);
    }

    for (const auto& ev_constraint : constraint_set.ev_constraints) {
        // ∇²[m * b(d(x))] = m * ∇(b'(d(x)) * ∇d(x))
        //                 = m * [b"(d(x)) * ∇d(x) * ∇d(x)ᵀ + b'(d(x)) * ∇²d(x)]
        const auto& p = V.row(ev_constraint.vertex_index);
        const auto& e0 = V.row(E(ev_constraint.edge_index, 0));
        const auto& e1 = V.row(E(ev_constraint.edge_index, 1));

        double distance_sqr =
            point_edge_distance(p, e0, e1, PointEdgeDistanceType::P_E);
        Eigen::VectorXd local_grad;
        point_edge_distance_gradient(
            p, e0, e1, PointEdgeDistanceType::P_E, local_grad);
        Eigen::MatrixXd local_hess;
        point_edge_distance_hessian(
            p, e0, e1, PointEdgeDistanceType::P_E, local_hess);

        local_hess *= barrier_gradient(distance_sqr, dhat_squared);
        local_hess += barrier_hessian(distance_sqr, dhat_squared) * local_grad
            * local_grad.transpose();

        local_hess *= ev_constraint.multiplicity;

        if (project_to_psd) {
            local_hess = Eigen::project_to_psd(local_hess);
        }

        // Map from local to global gradient
        std::vector<long> ids = { { ev_constraint.vertex_index,
                                    E(ev_constraint.edge_index, 0),
                                    E(ev_constraint.edge_index, 1) } };
        local_hessian_to_global_triplets(local_hess, ids, dim, hess_triplets);
    }

    for (const auto& ee_constraint : constraint_set.ee_constraints) {
        // ∇²[m(x) * b(d(x))] = ∇[∇m(x) * b(d(x)) + m(x) * b'(d(x)) * ∇d(x)]
        //                    = ∇²m(x) * b(d(x)) + b'(d(x)) * ∇d(x) * ∇m(x)ᵀ
        //                      + ∇m(x) * b'(d(x)) * ∇d(x))ᵀ
        //                      + m(x) * b"(d(x)) * ∇d(x) * ∇d(x)ᵀ
        //                      + m(x) * b'(d(x)) * ∇²d(x)
        const auto& ea0 = V.row(E(ee_constraint.edge0_index, 0));
        const auto& ea1 = V.row(E(ee_constraint.edge0_index, 1));
        const auto& eb0 = V.row(E(ee_constraint.edge1_index, 0));
        const auto& eb1 = V.row(E(ee_constraint.edge1_index, 1));

        // Compute distance derivatives
        // The distance type is unknown because of mollified PP and PE
        // constraints where also added as EE constraints.
        EdgeEdgeDistanceType dtype =
            edge_edge_distance_type(ea0, ea1, eb0, eb1);
        double distance_sqr = edge_edge_distance(ea0, ea1, eb0, eb1, dtype);
        Eigen::VectorXd distance_grad;
        edge_edge_distance_gradient(ea0, ea1, eb0, eb1, dtype, distance_grad);
        Eigen::MatrixXd distance_hess;
        edge_edge_distance_hessian(ea0, ea1, eb0, eb1, dtype, distance_hess);

        // Compute mollifier derivatives
        double mollifier =
            edge_edge_mollifier(ea0, ea1, eb0, eb1, ee_constraint.eps_x);
        Eigen::VectorXd mollifier_grad;
        edge_edge_mollifier_gradient(
            ea0, ea1, eb0, eb1, ee_constraint.eps_x, mollifier_grad);
        Eigen::MatrixXd mollifier_hess;
        edge_edge_mollifier_hessian(
            ea0, ea1, eb0, eb1, ee_constraint.eps_x, mollifier_hess);

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

        if (project_to_psd) {
            local_hess = Eigen::project_to_psd(local_hess);
        }

        // Map from local to global gradient
        std::vector<long> ids = {
            { E(ee_constraint.edge0_index, 0), E(ee_constraint.edge0_index, 1),
              E(ee_constraint.edge1_index, 0), E(ee_constraint.edge1_index, 1) }
        };
        local_hessian_to_global_triplets(local_hess, ids, dim, hess_triplets);
    }

    for (const auto& fv_constraint : constraint_set.fv_constraints) {
        // ∇²b(d(x)) = ∇(b'(d(x)) * ∇d(x))
        //           = b"(d(x)) * ∇d(x) * ∇d(x)ᵀ + b'(d(x)) * ∇²d(x)
        const auto& p = V.row(fv_constraint.vertex_index);
        const auto& t0 = V.row(F(fv_constraint.face_index, 0));
        const auto& t1 = V.row(F(fv_constraint.face_index, 1));
        const auto& t2 = V.row(F(fv_constraint.face_index, 2));

        double distance_sqr = point_triangle_distance(
            p, t0, t1, t2, PointTriangleDistanceType::P_T);
        Eigen::VectorXd local_grad;
        point_triangle_distance_gradient(
            p, t0, t1, t2, PointTriangleDistanceType::P_T, local_grad);
        Eigen::MatrixXd local_hess;
        point_triangle_distance_hessian(
            p, t0, t1, t2, PointTriangleDistanceType::P_T, local_hess);

        local_hess *= barrier_gradient(distance_sqr, dhat_squared);
        local_hess += barrier_hessian(distance_sqr, dhat_squared) * local_grad
            * local_grad.transpose();

        if (project_to_psd) {
            local_hess = Eigen::project_to_psd(local_hess);
        }

        // Map from local to global gradient
        std::vector<long> ids = {
            { fv_constraint.vertex_index, F(fv_constraint.face_index, 0),
              F(fv_constraint.face_index, 1), F(fv_constraint.face_index, 2) }
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
    bool ignore_internal_vertices,
    const Eigen::VectorXi& vertex_group_ids)
{
    int dim = V0.cols();
    assert(V1.cols() == dim);

    // Broad phase
    Candidates candidates;
    HashGrid hash_grid;
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
            E, vertex_group_ids, candidates.ev_candidates);
    } else {
        // These are not needed for 2D
        hash_grid.getEdgeEdgePairs(
            E, vertex_group_ids, candidates.ee_candidates);
        hash_grid.getFaceVertexPairs(
            F, vertex_group_ids, candidates.fv_candidates);
    }

    // Narrow phase
    double eta = 1e-6;

    for (const auto& ev_candidate : candidates.ev_candidates) {
        double toi;
        bool is_collision = point_edge_ccd(
            // Point at t=0
            V0.row(ev_candidate.vertex_index),
            // Edge at t=0
            V0.row(E(ev_candidate.edge_index, 0)),
            V0.row(E(ev_candidate.edge_index, 1)),
            // Point at t=1
            V1.row(ev_candidate.vertex_index),
            // Edge at t=1
            V1.row(E(ev_candidate.edge_index, 0)),
            V1.row(E(ev_candidate.edge_index, 1)), //
            toi);

        if (is_collision) {
            return false;
        }
    }

    for (const auto& ee_candidate : candidates.ee_candidates) {
        double toi;
        bool is_collision = edge_edge_ccd(
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
            toi);

        if (is_collision) {
            return false;
        }
    }

    for (const auto& fv_candidate : candidates.fv_candidates) {
        double toi;
        bool is_collision = point_triangle_ccd(
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
            toi);

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
    bool ignore_internal_vertices,
    const Eigen::VectorXi& vertex_group_ids)
{
    int dim = V0.cols();
    assert(V1.cols() == dim);

    // Broad phase
    Candidates candidates;
    HashGrid hash_grid;
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
            E, vertex_group_ids, candidates.ev_candidates);
    } else {
        // These are not needed for 2D
        hash_grid.getEdgeEdgePairs(
            E, vertex_group_ids, candidates.ee_candidates);
        hash_grid.getFaceVertexPairs(
            F, vertex_group_ids, candidates.fv_candidates);
    }

    // Narrow phase
    double earliest_toi = std::numeric_limits<double>::infinity();

    for (const auto& ev_candidate : candidates.ev_candidates) {
        double toi;
        bool is_collision = point_edge_ccd(
            // Point at t=0
            V0.row(ev_candidate.vertex_index),
            // Edge at t=0
            V0.row(E(ev_candidate.edge_index, 0)),
            V0.row(E(ev_candidate.edge_index, 1)),
            // Point at t=1
            V1.row(ev_candidate.vertex_index),
            // Edge at t=1
            V1.row(E(ev_candidate.edge_index, 0)),
            V1.row(E(ev_candidate.edge_index, 1)), //
            toi);

        assert(!is_collision || (toi >= 0 && toi <= 1));

        if (is_collision && toi < earliest_toi) {
            earliest_toi = toi;
        }
    }

    for (const auto& ee_candidate : candidates.ee_candidates) {
        double toi;
        bool is_collision = edge_edge_ccd(
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
            toi);

        assert(!is_collision || (toi >= 0 && toi <= 1));
        if (is_collision && toi < earliest_toi) {
            earliest_toi = toi;
        }
    }

    for (const auto& fv_candidate : candidates.fv_candidates) {
        double toi;
        bool is_collision = point_triangle_ccd(
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
            toi);

        assert(!is_collision || (toi >= 0 && toi <= 1));
        if (is_collision && toi < earliest_toi) {
            earliest_toi = toi;
        }
    }
    assert(earliest_toi >= 0);

    return std::min(earliest_toi, 1.0); // Fulfill the promise of a step size
}

// NOTE: Actually distance squared
double compute_minimum_distance(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const Constraints& constraint_set)
{
    double min_distance = std::numeric_limits<double>::infinity();

    for (const auto& vv_constraint : constraint_set.vv_constraints) {
        double distance_sqr = point_point_distance(
            V.row(vv_constraint.vertex0_index),
            V.row(vv_constraint.vertex1_index));

        if (distance_sqr < min_distance) {
            min_distance = distance_sqr;
        }
    }

    for (const auto& ev_constraint : constraint_set.ev_constraints) {
        double distance_sqr = point_edge_distance(
            V.row(ev_constraint.vertex_index),
            V.row(E(ev_constraint.edge_index, 0)),
            V.row(E(ev_constraint.edge_index, 1)), PointEdgeDistanceType::P_E);

        if (distance_sqr < min_distance) {
            min_distance = distance_sqr;
        }
    }

    for (const auto& ee_constraint : constraint_set.ee_constraints) {
        // The distance type is unknown because of mollified PP and PE
        // constraints where also added as EE constraints.
        double distance_sqr = edge_edge_distance(
            V.row(E(ee_constraint.edge0_index, 0)),
            V.row(E(ee_constraint.edge0_index, 1)),
            V.row(E(ee_constraint.edge1_index, 0)),
            V.row(E(ee_constraint.edge1_index, 1)));

        if (distance_sqr < min_distance) {
            min_distance = distance_sqr;
        }
    }

    for (const auto& fv_constraint : constraint_set.fv_constraints) {
        double distance_sqr = point_triangle_distance(
            V.row(fv_constraint.vertex_index),
            V.row(F(fv_constraint.face_index, 0)),
            V.row(F(fv_constraint.face_index, 1)),
            V.row(F(fv_constraint.face_index, 2)),
            PointTriangleDistanceType::P_T);

        if (distance_sqr < min_distance) {
            min_distance = distance_sqr;
        }
    }

    return min_distance;
}

} // namespace ipc
