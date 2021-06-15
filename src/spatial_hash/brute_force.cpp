#include <ipc/spatial_hash/brute_force.hpp>

#include <algorithm> // std::min/max

#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range2d.h>

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

    if (detect_edge_vertex) {
        detect_edge_vertex_collision_candidates_brute_force(
            V0, V1, E, candidates.ev_candidates, perform_aabb_check,
            aabb_inflation_radius, group_ids);
    }

    if (detect_edge_edge) {
        detect_edge_edge_collision_candidates_brute_force(
            V0, V1, E, candidates.ee_candidates, perform_aabb_check,
            aabb_inflation_radius, group_ids);
    }

    if (detect_face_vertex) {
        detect_face_vertex_collision_candidates_brute_force(
            V0, V1, F, candidates.fv_candidates, perform_aabb_check,
            aabb_inflation_radius, group_ids);
    }
}

template <typename Candidate>
void merge_local_candidates(
    const tbb::enumerable_thread_specific<std::vector<Candidate>>& storages,
    std::vector<Candidate>& candidates)
{
    // size up the candidates
    size_t num_candidates = candidates.size();
    for (const auto& local_candidates : storages) {
        num_candidates += local_candidates.size();
    }
    // serial merge!
    candidates.reserve(num_candidates);
    for (const auto& local_candidates : storages) {
        candidates.insert(
            candidates.end(), local_candidates.begin(), local_candidates.end());
    }
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

    tbb::enumerable_thread_specific<std::vector<EdgeVertexCandidate>> storages;

    const bool check_group = group_ids.size() > 0;
    tbb::parallel_for(
        tbb::blocked_range2d<long>(0l, long(E.rows()), 0l, long(V0.rows())),
        [&](const tbb::blocked_range2d<long>& r) {
            auto& local_candidates = storages.local();

            // Loop over edges
            for (long ei = r.rows().begin(); ei < r.rows().end(); ei++) {
                // Loop over vertices
                for (long vi = r.cols().begin(); vi < r.cols().end(); vi++) {

                    // Check that the vertex is not an endpoint of the edge
                    bool is_endpoint = vi == E(ei, 0) || vi == E(ei, 1);
                    bool same_group = check_group
                        && (group_ids(vi) == group_ids(E(ei, 0))
                            || group_ids(vi) == group_ids(E(ei, 1)));

                    if (!is_endpoint && !same_group) {
                        bool aabb_intersect = !perform_aabb_check
                            || point_edge_aabb_ccd(
                                   V0.row(vi), V0.row(E(ei, 0)),
                                   V0.row(E(ei, 1)), V1.row(vi),
                                   V1.row(E(ei, 0)), V1.row(E(ei, 1)),
                                   aabb_inflation_radius);
                        if (aabb_intersect) {
                            local_candidates.emplace_back(ei, vi);
                        }
                    }
                }
            }
        });

    merge_local_candidates(storages, ev_candidates);
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

    tbb::enumerable_thread_specific<std::vector<EdgeEdgeCandidate>> storages;

    const bool check_group = group_ids.size() > 0;
    long num_edges = E.rows();
    tbb::parallel_for(
        tbb::blocked_range2d<long>(0l, num_edges, 0l, num_edges),
        [&](const tbb::blocked_range2d<long>& r) {
            auto& local_candidates = storages.local();
            // eai < ebi
            long eai_end = std::min(r.rows().end(), r.cols().end());

            for (long eai = r.rows().begin(); eai < eai_end; eai++) {
                // i < r.cols().end() → i + 1 <= r.cols().end()
                int ebi_begin = std::max(eai + 1, r.cols().begin());
                assert(ebi_begin > eai);

                for (long ebi = ebi_begin; ebi < r.cols().end(); ebi++) {
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
                            local_candidates.emplace_back(eai, ebi);
                        }
                    }
                }
            }
        });

    merge_local_candidates(storages, ee_candidates);
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

    tbb::enumerable_thread_specific<std::vector<FaceVertexCandidate>> storages;

    const bool check_group = group_ids.size() > 0;
    tbb::parallel_for(
        tbb::blocked_range2d<long>(0l, long(F.rows()), 0l, long(V0.rows())),
        [&](const tbb::blocked_range2d<long>& r) {
            auto& local_candidates = storages.local();

            // Loop over faces
            for (int fi = r.rows().begin(); fi < r.rows().end(); fi++) {
                // Loop over vertices
                for (int vi = r.cols().begin(); vi < r.cols().end(); vi++) {
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
                                   V0.row(F(fi, 0)), V0.row(F(fi, 1)),
                                   V0.row(F(fi, 2)),
                                   V1.row(vi), //
                                   V1.row(F(fi, 0)), V1.row(F(fi, 1)),
                                   V1.row(F(fi, 2)), aabb_inflation_radius);
                        if (aabb_intersect) {
                            local_candidates.emplace_back(fi, vi);
                        }
                    }
                }
            }
        });

    merge_local_candidates(storages, fv_candidates);
}

} // namespace ipc
