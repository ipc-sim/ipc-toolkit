#include <ipc/ipc.hpp>

#include <unordered_map>

#include <ipc/ccd/ccd.hpp>
#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/edge_edge_mollifier.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/distance/point_triangle.hpp>
#include <ipc/spatial_hash/hash_grid.hpp>
#include <ipc/utils/local_to_global.hpp>

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
    bool ignore_codimensional_vertices,
    const Eigen::VectorXi& vertex_group_ids,
    const Eigen::MatrixXi& F2E,
    double dmin)
{
    double inflation_radius = (dhat + dmin) / 2;

    Candidates candidates;
    HashGrid hash_grid;
    hash_grid.resize(V, V, E, inflation_radius);

    // Assumes the edges connect to all boundary vertices
    if (ignore_codimensional_vertices) {
        for (int e = 0; e < E.rows(); ++e) {
            const int e0 = E(e, 0);
            const int e1 = E(e, 1);
            hash_grid.addVertex(
                V.row(e0), V.row(e0), e0,
                /*inflation_radius=*/dhat + dmin);
            hash_grid.addVertex(
                V.row(e1), V.row(e1), e1,
                /*inflation_radius=*/dhat + dmin);
        }
    } else {
        hash_grid.addVertices(V, V, inflation_radius);
    }

    hash_grid.addEdges(V, V, E, inflation_radius);
    if (V.cols() == 3) {
        // These are not needed for 2D
        hash_grid.addFaces(V, V, F, inflation_radius);
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

    construct_constraint_set(
        candidates, V_rest, V, E, F, dhat, constraint_set, F2E, dmin);
}

void construct_constraint_set(
    const Candidates& candidates,
    const Eigen::MatrixXd& V_rest,
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    double dhat,
    Constraints& constraint_set,
    const Eigen::MatrixXi& F2E,
    double dmin)
{
    double dhat_squared = dhat * dhat;
    double dmin_squared = dmin * dmin;

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

        if (distance_sqr - dmin_squared < 2 * dmin * dhat + dhat_squared) {
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
                constraint_set.ev_constraints.emplace_back(ev_candidate);
                ev_to_index.emplace(
                    constraint_set.ev_constraints.back(),
                    constraint_set.ev_constraints.size() - 1);
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

        if (distance_sqr - dmin_squared < 2 * dmin * dhat + dhat_squared) {
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
        long vi = fv_candidate.vertex_index;
        long fi = fv_candidate.face_index;
        long f0i = F(fi, 0);
        long f1i = F(fi, 1);
        long f2i = F(fi, 2);

        // Compute distance type
        PointTriangleDistanceType dtype = point_triangle_distance_type(
            V.row(fv_candidate.vertex_index), //
            V.row(f0i), V.row(f1i), V.row(f2i));

        double distance_sqr = point_triangle_distance(
            V.row(fv_candidate.vertex_index), //
            V.row(f0i), V.row(f1i), V.row(f2i), dtype);

        if (distance_sqr - dmin_squared < 2 * dmin * dhat + dhat_squared) {
            switch (dtype) {
            case PointTriangleDistanceType::P_T0:
                add_vertex_vertex_constraint(
                    constraint_set.vv_constraints, vv_to_index, vi, f0i);
                break;

            case PointTriangleDistanceType::P_T1:
                add_vertex_vertex_constraint(
                    constraint_set.vv_constraints, vv_to_index, vi, f1i);
                break;

            case PointTriangleDistanceType::P_T2:
                add_vertex_vertex_constraint(
                    constraint_set.vv_constraints, vv_to_index, vi, f2i);
                break;

            case PointTriangleDistanceType::P_E0:
                add_edge_vertex_constraint(
                    constraint_set.ev_constraints, ev_to_index,
                    F2E.rows() > fi ? F2E(fi, 0) : find_edge(E, f0i, f1i), vi);
                break;

            case PointTriangleDistanceType::P_E1:
                add_edge_vertex_constraint(
                    constraint_set.ev_constraints, ev_to_index,
                    F2E.rows() > fi ? F2E(fi, 1) : find_edge(E, f1i, f2i), vi);
                break;

            case PointTriangleDistanceType::P_E2:
                add_edge_vertex_constraint(
                    constraint_set.ev_constraints, ev_to_index,
                    F2E.rows() > fi ? F2E(fi, 2) : find_edge(E, f2i, f0i), vi);
                break;

            case PointTriangleDistanceType::P_T:
                constraint_set.fv_constraints.emplace_back(fv_candidate);
                break;
            }
        }
    }

    for (int ci = 0; ci < constraint_set.size(); ci++) {
        constraint_set[ci].minimum_distance = dmin;
    }
}

///////////////////////////////////////////////////////////////////////////////

double compute_barrier_potential(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const Constraints& constraint_set,
    double dhat)
{
    double potential = 0;
    for (size_t i = 0; i < constraint_set.size(); i++) {
        potential += constraint_set[i].compute_potential(V, E, F, dhat);
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
    Eigen::VectorXd grad = Eigen::VectorXd::Zero(V.size());
    int dim = V.cols();

    for (size_t i = 0; i < constraint_set.size(); i++) {
        local_gradient_to_global_gradient(
            constraint_set[i].compute_potential_gradient(V, E, F, dhat),
            constraint_set[i].vertex_indices(E, F), dim, grad);
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

    int dim = V.cols();
    int dim_sq = dim * dim;
    std::vector<Eigen::Triplet<double>> hess_triplets;
    hess_triplets.reserve(
        constraint_set.vv_constraints.size() * /*2*2=*/4 * dim_sq
        + constraint_set.ev_constraints.size() * /*3*3=*/9 * dim_sq
        + constraint_set.ee_constraints.size() * /*4*4=*/16 * dim_sq
        + constraint_set.fv_constraints.size() * /*4*4=*/16 * dim_sq);

    for (size_t i = 0; i < constraint_set.size(); i++) {
        local_hessian_to_global_triplets(
            constraint_set[i].compute_potential_hessian(
                V, E, F, dhat, project_to_psd),
            constraint_set[i].vertex_indices(E, F), dim, hess_triplets);
    }

    Eigen::SparseMatrix<double> hess(V.size(), V.size());
    hess.setFromTriplets(hess_triplets.begin(), hess_triplets.end());
    return hess;
}

///////////////////////////////////////////////////////////////////////////////

bool is_step_collision_free(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    bool ignore_codimensional_vertices,
    const Eigen::VectorXi& vertex_group_ids)
{
    int dim = V0.cols();
    assert(V1.cols() == dim);

    // Broad phase
    Candidates candidates;
    HashGrid hash_grid;
    hash_grid.resize(V0, V1, E);

    // Assumes the edges connect to all boundary vertices
    if (ignore_codimensional_vertices) {
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

///////////////////////////////////////////////////////////////////////////////

double compute_collision_free_stepsize(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    bool ignore_codimensional_vertices,
    const Eigen::VectorXi& vertex_group_ids)
{
    int dim = V0.cols();
    assert(V1.cols() == dim);

    // Broad phase
    Candidates candidates;
    HashGrid hash_grid;
    hash_grid.resize(V0, V1, E);

    // Assumes the edges connect to all boundary vertices
    if (ignore_codimensional_vertices) {
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
    double earliest_toi = 1;

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
            toi,
            /*tmax=*/earliest_toi);

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
            toi,
            /*tmax=*/earliest_toi);

        assert(!is_collision || (toi >= 0 && toi <= 1));
        if (is_collision && toi < earliest_toi) {
            earliest_toi = toi;
        }
    }

    assert(earliest_toi >= 0 && earliest_toi <= 1.0);
    return earliest_toi;
}

///////////////////////////////////////////////////////////////////////////////

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
