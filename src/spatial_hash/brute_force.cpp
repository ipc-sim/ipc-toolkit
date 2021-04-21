#include <ipc/spatial_hash/brute_force.hpp>

#include <tbb/parallel_invoke.h>

#include <ipc/ccd/broadphase.hpp>

namespace ipc {

// Find all edge-vertex collisions in one time step using brute-force
// comparisons of all E and all vertices.
void detect_collision_candidates_brute_force(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    Candidates& candidates,
    bool detect_edge_vertex,
    bool detect_edge_edge,
    bool detect_face_vertex,
    bool perform_aabb_check,
    double aabb_inflation_radius,
    const Eigen::VectorXi& group_ids)
{
    assert(V0.rows() == V1.rows() && V0.cols() == V1.cols());
    assert(E.size() == 0 || E.cols() == 2);
    assert(F.size() == 0 || F.cols() == 3);

    // Loop over all E
    tbb::parallel_invoke(
        [&] {
            if (detect_edge_vertex) {
                detect_edge_vertex_collision_candidates_brute_force(
                    V0, V1, E, candidates.ev_candidates, perform_aabb_check,
                    aabb_inflation_radius, group_ids);
            }
        },
        [&] {
            if (detect_edge_edge) {
                detect_edge_edge_collision_candidates_brute_force(
                    V0, V1, E, candidates.ee_candidates, perform_aabb_check,
                    aabb_inflation_radius, group_ids);
            }
        },
        [&] {
            if (detect_face_vertex) {
                detect_face_vertex_collision_candidates_brute_force(
                    V0, V1, F, candidates.fv_candidates, perform_aabb_check,
                    aabb_inflation_radius, group_ids);
            }
        });
}

void detect_edge_vertex_collision_candidates_brute_force(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    std::vector<EdgeVertexCandidate>& ev_candidates,
    bool perform_aabb_check,
    double aabb_inflation_radius,
    const Eigen::VectorXi& group_ids)
{
    assert(V0.rows() == V1.rows() && V0.cols() == V1.cols());

    const bool check_group = group_ids.size() > 0;
    for (int ei = 0; ei < E.rows(); ei++) {
        // Loop over all vertices
        for (int vi = 0; vi < V0.rows(); vi++) {
            // Check that the vertex is not an endpoint of the edge
            bool is_endpoint = vi == E(ei, 0) || vi == E(ei, 1);
            bool same_group = check_group
                && (group_ids(vi) == group_ids(E(ei, 0))
                    || group_ids(vi) == group_ids(E(ei, 1)));
            if (!is_endpoint && !same_group) {
                bool aabb_intersect = !perform_aabb_check
                    || point_edge_aabb_ccd(
                           V0.row(vi), V0.row(E(ei, 0)), V0.row(E(ei, 1)),
                           V1.row(vi), V1.row(E(ei, 0)), V1.row(E(ei, 1)),
                           aabb_inflation_radius);
                if (aabb_intersect) {
                    ev_candidates.emplace_back(ei, vi);
                }
            }
        }
    }
}

void detect_edge_edge_collision_candidates_brute_force(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    std::vector<EdgeEdgeCandidate>& ee_candidates,
    bool perform_aabb_check,
    double aabb_inflation_radius,
    const Eigen::VectorXi& group_ids)
{
    assert(V0.rows() == V1.rows() && V0.cols() == V1.cols());

    const bool check_group = group_ids.size() > 0;
    for (int eai = 0; eai < E.rows(); eai++) {
        // Loop over all remaining E
        for (int ebi = eai + 1; ebi < E.rows(); ebi++) {
            bool has_common_endpoint = E(eai, 0) == E(ebi, 0)
                || E(eai, 0) == E(ebi, 1) || E(eai, 1) == E(ebi, 0)
                || E(eai, 1) == E(ebi, 1);
            bool same_group = check_group
                && (group_ids(E(eai, 0)) == group_ids(E(ebi, 0))
                    || group_ids(E(eai, 0)) == group_ids(E(ebi, 1))
                    || group_ids(E(eai, 1)) == group_ids(E(ebi, 0))
                    || group_ids(E(eai, 1)) == group_ids(E(ebi, 1)));
            if (!has_common_endpoint && !same_group) {
                bool aabb_intersect = !perform_aabb_check
                    || edge_edge_aabb_ccd(
                           V0.row(E(eai, 0)), V0.row(E(eai, 1)), //
                           V0.row(E(ebi, 0)), V0.row(E(ebi, 1)), //
                           V1.row(E(eai, 0)), V1.row(E(eai, 1)), //
                           V1.row(E(ebi, 0)), V1.row(E(ebi, 1)), //
                           aabb_inflation_radius);
                if (aabb_intersect) {
                    ee_candidates.emplace_back(eai, ebi);
                }
            }
        }
    }
}

void detect_face_vertex_collision_candidates_brute_force(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& F,
    std::vector<FaceVertexCandidate>& fv_candidates,
    bool perform_aabb_check,
    double aabb_inflation_radius,
    const Eigen::VectorXi& group_ids)
{
    assert(V0.rows() == V1.rows() && V0.cols() == V1.cols());

    const bool check_group = group_ids.size() > 0;
    // Loop over all F
    for (int fi = 0; fi < F.rows(); fi++) {
        // Loop over all vertices
        for (int vi = 0; vi < V0.rows(); vi++) {
            // Check that the vertex is not an endpoint of the edge
            bool is_endpoint =
                vi == F(fi, 0) || vi == F(fi, 1) || vi == F(fi, 2);
            bool same_group = check_group
                && (group_ids(vi) == group_ids(F(fi, 0))
                    || group_ids(vi) == group_ids(F(fi, 1))
                    || group_ids(vi) == group_ids(F(fi, 2)));
            if (!is_endpoint && !same_group) {
                bool aabb_intersect = !perform_aabb_check
                    || point_triangle_aabb_ccd(
                           V0.row(vi), //
                           V0.row(F(fi, 0)), V0.row(F(fi, 1)), V0.row(F(fi, 2)),
                           V1.row(vi), //
                           V1.row(F(fi, 0)), V1.row(F(fi, 1)), V1.row(F(fi, 2)),
                           aabb_inflation_radius);
                if (aabb_intersect) {
                    fv_candidates.emplace_back(fi, vi);
                }
            }
        }
    }
}

} // namespace ipc
