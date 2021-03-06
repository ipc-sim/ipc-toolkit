#include <ipc/broad_phase/brute_force.hpp>

#include <algorithm> // std::min/max

#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range2d.h>

#include <ipc/ccd/aabb.hpp>

namespace ipc {

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

///////////////////////////////////////////////////////////////////////////////
// Discrete Collision Detection

void build_is_codim(
    size_t num_vertices, const Eigen::MatrixXi& E, std::vector<bool>& is_codim)
{
    is_codim.resize(num_vertices, true);
    // Column first because colmajor
    for (size_t ej = 0; ej < E.cols(); ej++) {
        for (size_t ei = 0; ei < E.rows(); ei++) {
            assert(E(ei, ej) < num_vertices);
            is_codim[E(ei, ej)] = false;
        }
    }
}

void detect_collision_candidates_brute_force(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    Candidates& candidates,
    bool detect_edge_vertex,
    bool detect_edge_edge,
    bool detect_face_vertex,
    bool perform_aabb_check,
    double aabb_inflation_radius,
    bool ignore_codimensional_vertices,
    const std::function<bool(size_t, size_t)>& can_collide)
{
    assert(E.size() == 0 || E.cols() == 2);
    assert(F.size() == 0 || F.cols() == 3);

    std::vector<bool> is_codim;
    if (ignore_codimensional_vertices) {
        build_is_codim(V.rows(), E, is_codim);
    }
    auto adjusted_can_collide = [&](size_t vi, size_t vj) {
        return (is_codim.empty() || (!is_codim[vi] && !is_codim[vj]))
            && can_collide(vi, vj);
    };

    if (detect_edge_vertex) {
        detect_edge_vertex_collision_candidates_brute_force(
            V, E, candidates.ev_candidates, perform_aabb_check,
            aabb_inflation_radius, adjusted_can_collide);
    }

    if (detect_edge_edge) {
        // Use the original can_collide because edge vertices are not
        // codimensional.
        detect_edge_edge_collision_candidates_brute_force(
            V, E, candidates.ee_candidates, perform_aabb_check,
            aabb_inflation_radius, can_collide);
    }

    if (detect_face_vertex) {
        detect_face_vertex_collision_candidates_brute_force(
            V, F, candidates.fv_candidates, perform_aabb_check,
            aabb_inflation_radius, adjusted_can_collide);
    }
}

void detect_edge_vertex_collision_candidates_brute_force(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    std::vector<EdgeVertexCandidate>& ev_candidates,
    bool perform_aabb_check,
    double aabb_inflation_radius,
    const std::function<bool(size_t, size_t)>& can_collide)
{
    tbb::enumerable_thread_specific<std::vector<EdgeVertexCandidate>> storages;

    tbb::parallel_for(
        tbb::blocked_range2d<long>(0l, long(E.rows()), 0l, long(V.rows())),
        [&](const tbb::blocked_range2d<long>& r) {
            auto& local_candidates = storages.local();

            // Loop over edges
            for (long ei = r.rows().begin(); ei < r.rows().end(); ei++) {
                const size_t &e0i = E(ei, 0), &e1i = E(ei, 1);

                // Loop over vertices
                for (long vi = r.cols().begin(); vi < r.cols().end(); vi++) {
                    // Check that the vertex is not an endpoint of the edge
                    if (vi == e0i || vi == e1i) {
                        continue;
                    }

                    if (!can_collide(vi, e0i) && !can_collide(vi, e1i)) {
                        continue;
                    }

                    bool aabb_intersect = !perform_aabb_check
                        || point_edge_aabb_cd(
                               V.row(vi), V.row(E(ei, 0)), V.row(E(ei, 1)),
                               2 * aabb_inflation_radius);
                    if (aabb_intersect) {
                        local_candidates.emplace_back(ei, vi);
                    }
                }
            }
        });

    merge_local_candidates(storages, ev_candidates);
}

void detect_edge_edge_collision_candidates_brute_force(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    std::vector<EdgeEdgeCandidate>& ee_candidates,
    bool perform_aabb_check,
    double aabb_inflation_radius,
    const std::function<bool(size_t, size_t)>& can_collide)
{
    tbb::enumerable_thread_specific<std::vector<EdgeEdgeCandidate>> storages;

    long num_edges = E.rows();
    tbb::parallel_for(
        tbb::blocked_range2d<long>(0l, num_edges, 0l, num_edges),
        [&](const tbb::blocked_range2d<long>& r) {
            auto& local_candidates = storages.local();
            // eai < ebi
            long eai_end = std::min(r.rows().end(), r.cols().end());

            for (long eai = r.rows().begin(); eai < eai_end; eai++) {
                const size_t &ea0i = E(eai, 0), &ea1i = E(eai, 1);

                // i < r.cols().end() → i + 1 <= r.cols().end()
                int ebi_begin = std::max(eai + 1, r.cols().begin());
                assert(ebi_begin > eai);

                for (long ebi = ebi_begin; ebi < r.cols().end(); ebi++) {
                    const size_t &eb0i = E(ebi, 0), &eb1i = E(ebi, 1);

                    // Check for common end points
                    if (ea0i == eb0i || ea0i == eb1i || ea1i == eb0i
                        || ea1i == eb1i) {
                        continue;
                    }

                    if (!can_collide(ea0i, eb0i) && !can_collide(ea0i, eb1i)
                        && !can_collide(ea1i, eb0i)
                        && !can_collide(ea1i, eb1i)) {
                        continue;
                    }

                    bool aabb_intersect = !perform_aabb_check
                        || edge_edge_aabb_cd(
                               V.row(ea0i), V.row(ea1i), V.row(eb0i),
                               V.row(eb1i), 2 * aabb_inflation_radius);
                    if (aabb_intersect) {
                        local_candidates.emplace_back(eai, ebi);
                    }
                }
            }
        });

    merge_local_candidates(storages, ee_candidates);
}

void detect_face_vertex_collision_candidates_brute_force(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    std::vector<FaceVertexCandidate>& fv_candidates,
    bool perform_aabb_check,
    double aabb_inflation_radius,
    const std::function<bool(size_t, size_t)>& can_collide)
{
    tbb::enumerable_thread_specific<std::vector<FaceVertexCandidate>> storages;

    tbb::parallel_for(
        tbb::blocked_range2d<long>(0l, long(F.rows()), 0l, long(V.rows())),
        [&](const tbb::blocked_range2d<long>& r) {
            auto& local_candidates = storages.local();

            // Loop over faces
            for (int fi = r.rows().begin(); fi < r.rows().end(); fi++) {
                const size_t &f0i = F(fi, 0), &f1i = F(fi, 1), &f2i = F(fi, 2);

                // Loop over vertices
                for (int vi = r.cols().begin(); vi < r.cols().end(); vi++) {
                    // Check for common end points
                    if (vi == f0i || vi == f1i || vi == f2i) {
                        continue;
                    }

                    if (!can_collide(vi, f0i) && !can_collide(vi, f1i)
                        && !can_collide(vi, f2i)) {
                        continue;
                    }

                    bool aabb_intersect = !perform_aabb_check
                        || point_triangle_aabb_cd(
                               V.row(vi), V.row(f0i), V.row(f1i), V.row(f2i),
                               2 * aabb_inflation_radius);
                    if (aabb_intersect) {
                        local_candidates.emplace_back(fi, vi);
                    }
                }
            }
        });

    merge_local_candidates(storages, fv_candidates);
}

void detect_edge_face_collision_candidates_brute_force(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    std::vector<EdgeFaceCandidate>& ef_candidates,
    bool perform_aabb_check,
    double aabb_inflation_radius,
    const std::function<bool(size_t, size_t)>& can_collide)
{
    tbb::enumerable_thread_specific<std::vector<EdgeFaceCandidate>> storages;

    tbb::parallel_for(
        tbb::blocked_range2d<long>(0l, long(E.rows()), 0l, long(F.rows())),
        [&](const tbb::blocked_range2d<long>& r) {
            auto& local_candidates = storages.local();

            // Loop over faces
            for (int ei = r.rows().begin(); ei < r.rows().end(); ei++) {
                const size_t &e0i = E(ei, 0), &e1i = E(ei, 1);

                // Loop over vertices
                for (int fi = r.cols().begin(); fi < r.cols().end(); fi++) {
                    const size_t &f0i = F(fi, 0), &f1i = F(fi, 1),
                                 &f2i = F(fi, 2);

                    // Check for common end points
                    if (e0i == f0i || e0i == f1i || e0i == f2i || e1i == f0i
                        || e1i == f1i || e1i == f2i) {
                        continue;
                    }

                    if (!can_collide(e0i, f0i) && !can_collide(e0i, f1i)
                        && !can_collide(e0i, f2i) && !can_collide(e1i, f0i)
                        && !can_collide(e1i, f1i) && !can_collide(e1i, f2i)) {
                        continue;
                    }

                    bool aabb_intersect = !perform_aabb_check
                        || edge_triangle_aabb_cd(
                               V.row(e0i), V.row(e1i), V.row(f0i), V.row(f1i),
                               V.row(f2i), 2 * aabb_inflation_radius);
                    if (aabb_intersect) {
                        local_candidates.emplace_back(ei, fi);
                    }
                }
            }
        });

    merge_local_candidates(storages, ef_candidates);
}

///////////////////////////////////////////////////////////////////////////////
// Continuous Collision Detection

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
    bool ignore_codimensional_vertices,
    const std::function<bool(size_t, size_t)>& can_collide)
{
    assert(V0.rows() == V1.rows() && V0.cols() == V1.cols());
    assert(E.size() == 0 || E.cols() == 2);
    assert(F.size() == 0 || F.cols() == 3);

    std::vector<bool> is_codim;
    if (ignore_codimensional_vertices) {
        build_is_codim(V0.rows(), E, is_codim);
    }
    auto adjusted_can_collide = [&](size_t vi, size_t vj) {
        return (is_codim.empty() || (!is_codim[vi] && !is_codim[vj]))
            && can_collide(vi, vj);
    };

    if (detect_edge_vertex) {
        detect_edge_vertex_collision_candidates_brute_force(
            V0, V1, E, candidates.ev_candidates, perform_aabb_check,
            aabb_inflation_radius, adjusted_can_collide);
    }

    if (detect_edge_edge) {
        // Use the original can_collide because edge vertices are not
        // codimensional.
        detect_edge_edge_collision_candidates_brute_force(
            V0, V1, E, candidates.ee_candidates, perform_aabb_check,
            aabb_inflation_radius, can_collide);
    }

    if (detect_face_vertex) {
        detect_face_vertex_collision_candidates_brute_force(
            V0, V1, F, candidates.fv_candidates, perform_aabb_check,
            aabb_inflation_radius, adjusted_can_collide);
    }
}

void detect_edge_vertex_collision_candidates_brute_force(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& E,
    std::vector<EdgeVertexCandidate>& ev_candidates,
    bool perform_aabb_check,
    double aabb_inflation_radius,
    const std::function<bool(size_t, size_t)>& can_collide)
{
    assert(V0.rows() == V1.rows() && V0.cols() == V1.cols());

    tbb::enumerable_thread_specific<std::vector<EdgeVertexCandidate>> storages;

    tbb::parallel_for(
        tbb::blocked_range2d<long>(0l, long(E.rows()), 0l, long(V0.rows())),
        [&](const tbb::blocked_range2d<long>& r) {
            auto& local_candidates = storages.local();

            // Loop over edges
            for (long ei = r.rows().begin(); ei < r.rows().end(); ei++) {
                const size_t &e0i = E(ei, 0), &e1i = E(ei, 1);

                // Loop over vertices
                for (long vi = r.cols().begin(); vi < r.cols().end(); vi++) {
                    // Check that the vertex is not an endpoint of the edge
                    if (vi == e0i || vi == e1i) {
                        continue;
                    }

                    if (!can_collide(vi, e0i) && !can_collide(vi, e1i)) {
                        continue;
                    }

                    bool aabb_intersect = !perform_aabb_check
                        || point_edge_aabb_ccd(
                               V0.row(vi), V0.row(e0i), V0.row(e0i), V1.row(vi),
                               V1.row(e0i), V1.row(e1i),
                               2 * aabb_inflation_radius);
                    if (aabb_intersect) {
                        local_candidates.emplace_back(ei, vi);
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
    const std::function<bool(size_t, size_t)>& can_collide)
{
    assert(V0.rows() == V1.rows() && V0.cols() == V1.cols());

    tbb::enumerable_thread_specific<std::vector<EdgeEdgeCandidate>> storages;

    long num_edges = E.rows();
    tbb::parallel_for(
        tbb::blocked_range2d<long>(0l, num_edges, 0l, num_edges),
        [&](const tbb::blocked_range2d<long>& r) {
            auto& local_candidates = storages.local();
            // eai < ebi
            long eai_end = std::min(r.rows().end(), r.cols().end());

            for (long eai = r.rows().begin(); eai < eai_end; eai++) {
                const size_t &ea0i = E(eai, 0), &ea1i = E(eai, 1);

                // i < r.cols().end() → i + 1 <= r.cols().end()
                int ebi_begin = std::max(eai + 1, r.cols().begin());
                assert(ebi_begin > eai);

                for (long ebi = ebi_begin; ebi < r.cols().end(); ebi++) {
                    const size_t &eb0i = E(ebi, 0), &eb1i = E(ebi, 1);

                    // Check for common end points
                    if (ea0i == eb0i || ea0i == eb1i || ea1i == eb0i
                        || ea1i == eb1i) {
                        continue;
                    }

                    if (!can_collide(ea0i, eb0i) && !can_collide(ea0i, eb1i)
                        && !can_collide(ea1i, eb0i)
                        && !can_collide(ea1i, eb1i)) {
                        continue;
                    }

                    bool aabb_intersect = !perform_aabb_check
                        || edge_edge_aabb_ccd(
                               V0.row(ea0i), V0.row(ea1i), V0.row(eb0i),
                               V0.row(eb1i), V1.row(ea0i), V1.row(ea1i),
                               V1.row(eb0i), V1.row(eb1i),
                               2 * aabb_inflation_radius);
                    if (aabb_intersect) {
                        local_candidates.emplace_back(eai, ebi);
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
    const std::function<bool(size_t, size_t)>& can_collide)
{
    assert(V0.rows() == V1.rows() && V0.cols() == V1.cols());

    tbb::enumerable_thread_specific<std::vector<FaceVertexCandidate>> storages;

    tbb::parallel_for(
        tbb::blocked_range2d<long>(0l, long(F.rows()), 0l, long(V0.rows())),
        [&](const tbb::blocked_range2d<long>& r) {
            auto& local_candidates = storages.local();

            // Loop over faces
            for (int fi = r.rows().begin(); fi < r.rows().end(); fi++) {
                const size_t &f0i = F(fi, 0), &f1i = F(fi, 1), &f2i = F(fi, 2);

                // Loop over vertices
                for (int vi = r.cols().begin(); vi < r.cols().end(); vi++) {
                    // Check for common end points
                    if (vi == f0i || vi == f1i || vi == f2i) {
                        continue;
                    }

                    if (!can_collide(vi, f0i) && !can_collide(vi, f1i)
                        && !can_collide(vi, f2i)) {
                        continue;
                    }

                    bool aabb_intersect = !perform_aabb_check
                        || point_triangle_aabb_ccd(
                               V0.row(vi), V0.row(f0i), V0.row(f1i),
                               V0.row(f2i), V1.row(vi), V1.row(f0i),
                               V1.row(f1i), V1.row(f2i),
                               2 * aabb_inflation_radius);
                    if (aabb_intersect) {
                        local_candidates.emplace_back(fi, vi);
                    }
                }
            }
        });

    merge_local_candidates(storages, fv_candidates);
}

} // namespace ipc
