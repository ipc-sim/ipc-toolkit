#include "candidates.hpp"

#include <ipc/config.hpp>
#include <ipc/ipc.hpp>
#include <ipc/broad_phase/default_broad_phase.hpp>
#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/distance/point_line.hpp>
#include <ipc/distance/point_plane.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/utils/eigen_ext.hpp>
#include <ipc/utils/profiler.hpp>
#include <ipc/utils/save_obj.hpp>
#include <ipc/utils/unordered_map_and_set.hpp>

#include <igl/remove_unreferenced.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/parallel_sort.h>

#include <atomic>
#include <fstream>
#include <iostream>

namespace ipc {

// Definition of the pimpl declared in candidates.hpp. Kept here so the public
// header does not need to include Abseil (a private dependency).
struct Candidates::AdjacencySets {
    unordered_map<index_t, std::set<index_t>> vv;
    unordered_map<index_t, std::set<index_t>> ve;
    unordered_map<index_t, std::set<index_t>> vf;

    unordered_map<index_t, std::set<index_t>> ev;
    unordered_map<index_t, std::set<index_t>> ee;
    unordered_map<index_t, std::set<index_t>> ef;

    unordered_map<index_t, std::set<index_t>> fv;
    unordered_map<index_t, std::set<index_t>> fe;
    unordered_map<index_t, std::set<index_t>> ff;
};

namespace {
    // Pad codim_edges because remove_unreferenced requires a N×3 matrix.
    Eigen::MatrixXi pad_edges(Eigen::ConstRef<Eigen::MatrixXi> E)
    {
        assert(E.cols() == 2);
        Eigen::MatrixXi E_padded(E.rows(), 3);
        E_padded.leftCols(2) = E;
        E_padded.col(2) = E.col(1);
        return E_padded;
    }

    Eigen::MatrixXi unpad_edges(Eigen::ConstRef<Eigen::MatrixXi> E_padded)
    {
        return E_padded.leftCols(2);
    }
} // namespace

void Candidates::build(
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices,
    const double inflation_radius,
    BroadPhase* broad_phase,
    const bool all_types)
{
    IPC_TOOLKIT_PROFILE_BLOCK("Candidates::build(static)");

    std::unique_ptr<BroadPhase> default_broad_phase;
    if (broad_phase == nullptr) {
        default_broad_phase = make_default_broad_phase();
        broad_phase = default_broad_phase.get();
    }

    const int dim = vertices.cols();
    m_mesh = mesh;

    clear();

    broad_phase->can_vertices_collide = mesh.can_collide;
    broad_phase->build(vertices, mesh.edges(), mesh.faces(), inflation_radius);
    broad_phase->detect_collision_candidates(*this, all_types);

    // Codim. vertices to codim. vertices:
    if (mesh.num_codim_vertices()) {
        broad_phase->clear();
        broad_phase->build(
            vertices(mesh.codim_vertices(), Eigen::all), //
            Eigen::MatrixXi(), Eigen::MatrixXi(), inflation_radius);

        broad_phase->detect_vertex_vertex_candidates(vv_candidates);
        for (auto& [vi, vj] : vv_candidates) {
            vi = mesh.codim_vertices()[vi];
            vj = mesh.codim_vertices()[vj];
        }
    }

    // Codim. edges to codim. vertices:
    // Only need this in 3D because in 2D, the codim. edges are the same as the
    // edges of the boundary. Only need codim. edge to codim. vertex because
    // codim. edge to non-codim. vertex is the same as edge-edge or face-vertex.
    if (dim == 3 && mesh.num_codim_vertices() && mesh.num_codim_edges()) {
        // Extract the vertices of the codim. edges
        Eigen::MatrixXd CE_V; // vertices of codim. edges
        Eigen::MatrixXi CE;   // codim. edges (indices into CEV)
        {
            Eigen::VectorXi _I, _J; // unused mappings
            igl::remove_unreferenced(
                vertices,
                pad_edges(mesh.edges()(mesh.codim_edges(), Eigen::all)), CE_V,
                CE, _I, _J);
            CE = unpad_edges(CE);
        }

        const size_t nCV = mesh.num_codim_vertices();
        Eigen::MatrixXd V(nCV + CE_V.rows(), dim);
        V.topRows(nCV) = vertices(mesh.codim_vertices(), Eigen::all);
        V.bottomRows(CE_V.rows()) = CE_V;

        CE.array() += nCV; // Offset indices to account for codim. vertices

        // TODO: Can we reuse the broad phase from above?
        broad_phase->clear();
        // Ignore c-edge to c-edge and c-vertex to c-vertex
        broad_phase->can_vertices_collide =
            make_codim_cross_filter(nCV) & mesh.can_collide;
        broad_phase->build(V, CE, Eigen::MatrixXi(), inflation_radius);

        broad_phase->detect_edge_vertex_candidates(ev_candidates);
        for (auto& [ei, vi] : ev_candidates) {
            assert(vi < mesh.codim_vertices().size());
            ei = mesh.codim_edges()[ei];    // Map back to mesh.edges
            vi = mesh.codim_vertices()[vi]; // Map back to vertices
        }
    }

    // Planes to vertices:
    for (const auto& plane : mesh.planes) {
        for (index_t vi = 0; vi < mesh.num_vertices(); ++vi) {
            if (plane.signedDistance(vertices.row(vi)) < inflation_radius) {
                pv_candidates.emplace_back(plane, vi);
            }
        }
    }
}

void Candidates::build(
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices_t0,
    Eigen::ConstRef<Eigen::MatrixXd> vertices_t1,
    const double inflation_radius,
    BroadPhase* broad_phase,
    const bool all_types)
{
    IPC_TOOLKIT_PROFILE_BLOCK("Candidates::build(dynamic)");

    std::unique_ptr<BroadPhase> default_broad_phase;
    if (broad_phase == nullptr) {
        default_broad_phase = make_default_broad_phase();
        broad_phase = default_broad_phase.get();
    }

    const int dim = vertices_t0.cols();
    m_mesh = mesh;

    clear();

    broad_phase->can_vertices_collide = mesh.can_collide;
    broad_phase->build(
        vertices_t0, vertices_t1, mesh.edges(), mesh.faces(), inflation_radius);
    broad_phase->detect_collision_candidates(*this);

    // Codim. vertices to codim. vertices:
    if (mesh.num_codim_vertices()) {
        broad_phase->clear();
        broad_phase->build(
            vertices_t0(mesh.codim_vertices(), Eigen::all),
            vertices_t1(mesh.codim_vertices(), Eigen::all), //
            Eigen::MatrixXi(), Eigen::MatrixXi(), inflation_radius);

        broad_phase->detect_vertex_vertex_candidates(vv_candidates);
        for (auto& [vi, vj] : vv_candidates) {
            vi = mesh.codim_vertices()[vi];
            vj = mesh.codim_vertices()[vj];
        }
    }

    // Codim. edges to codim. vertices:
    // Only need this in 3D because in 2D, the codim. edges are the same as the
    // edges of the boundary. Only need codim. edge to codim. vertex because
    // codim. edge to non-codim. vertex is the same as edge-edge or face-vertex.
    if (dim == 3 && mesh.num_codim_vertices() && mesh.num_codim_edges()) {
        // Extract the vertices of the codim. edges
        Eigen::MatrixXd CE_V_t0, CE_V_t1; // vertices of codim. edges
        Eigen::MatrixXi CE;               // codim. edges (indices into CEV)
        {
            Eigen::VectorXi _I, J;
            igl::remove_unreferenced(
                vertices_t0,
                pad_edges(mesh.edges()(mesh.codim_edges(), Eigen::all)),
                CE_V_t0, CE, _I, J);
            CE_V_t1 = vertices_t1(J, Eigen::all);
            CE = unpad_edges(CE);
        }

        const size_t nCV = mesh.num_codim_vertices();

        Eigen::MatrixXd V_t0(nCV + CE_V_t0.rows(), dim);
        V_t0.topRows(nCV) = vertices_t0(mesh.codim_vertices(), Eigen::all);
        V_t0.bottomRows(CE_V_t0.rows()) = CE_V_t0;

        Eigen::MatrixXd V_t1(nCV + CE_V_t1.rows(), dim);
        V_t1.topRows(nCV) = vertices_t1(mesh.codim_vertices(), Eigen::all);
        V_t1.bottomRows(CE_V_t1.rows()) = CE_V_t1;

        CE.array() += nCV; // Offset indices to account for codim. vertices

        // TODO: Can we reuse the broad phase from above?
        broad_phase->clear();
        // Ignore c-edge to c-edge and c-vertex to c-vertex
        broad_phase->can_vertices_collide =
            make_codim_cross_filter(nCV) & mesh.can_collide;
        broad_phase->build(V_t0, V_t1, CE, Eigen::MatrixXi(), inflation_radius);

        broad_phase->detect_edge_vertex_candidates(ev_candidates);
        for (auto& [ei, vi] : ev_candidates) {
            assert(vi < mesh.codim_vertices().size());
            ei = mesh.codim_edges()[ei];    // Map back to mesh.edges
            vi = mesh.codim_vertices()[vi]; // Map back to vertices
        }
    }

    // Planes to vertices:
    for (const auto& plane : mesh.planes) {
        for (index_t vi = 0; vi < mesh.num_vertices(); ++vi) {
            const double d_t0 = plane.signedDistance(vertices_t0.row(vi));
            const double d_t1 = plane.signedDistance(vertices_t1.row(vi));
            if (d_t0 < inflation_radius || d_t1 < inflation_radius) {
                pv_candidates.emplace_back(plane, vi);
            }
        }
    }
}

bool Candidates::is_step_collision_free(
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices_t0,
    Eigen::ConstRef<Eigen::MatrixXd> vertices_t1,
    const double min_distance,
    const NarrowPhaseCCD& narrow_phase_ccd) const
{
    assert(vertices_t0.rows() == mesh.num_vertices());
    assert(vertices_t1.rows() == mesh.num_vertices());

    // Narrow phase
    for (size_t i = 0; i < size(); i++) {
        const CollisionStencil& candidate = (*this)[i];

        double toi;
        bool is_collision = candidate.ccd(
            candidate.dof(vertices_t0, mesh.edges(), mesh.faces()),
            candidate.dof(vertices_t1, mesh.edges(), mesh.faces()), //
            toi, min_distance, /*tmax=*/1.0, narrow_phase_ccd);

        if (is_collision) {
            return false;
        }
    }

    return true;
}

double Candidates::compute_collision_free_stepsize(
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices_t0,
    Eigen::ConstRef<Eigen::MatrixXd> vertices_t1,
    const double min_distance,
    const NarrowPhaseCCD& narrow_phase_ccd) const
{
    assert(vertices_t0.rows() == mesh.num_vertices());
    assert(vertices_t1.rows() == mesh.num_vertices());
    IPC_TOOLKIT_PROFILE_BLOCK("Candidates::compute_collision_free_stepsize");

    if (empty()) {
        return 1; // No possible collisions, so can take full step.
    }

    std::atomic<double> earliest_toi(1.0);

    tbb::parallel_for(size_t(0), size(), [&](size_t i) {
        double tmax = earliest_toi.load(std::memory_order_relaxed);

        const CollisionStencil& candidate = (*this)[i];

        double toi = std::numeric_limits<double>::infinity(); // output
        const bool are_colliding = candidate.ccd(
            candidate.dof(vertices_t0, mesh.edges(), mesh.faces()),
            candidate.dof(vertices_t1, mesh.edges(), mesh.faces()), //
            toi, min_distance, tmax, narrow_phase_ccd);

        if (are_colliding) {
            // Update the earliest time of impact (TOI) atomically
            double prev = earliest_toi.load(std::memory_order_relaxed);
            while (toi < prev
                   && !earliest_toi.compare_exchange_weak(
                       prev, toi, std::memory_order_relaxed)) { }
        }
    });

    double result = earliest_toi.load(std::memory_order_relaxed);
    assert(result >= 0 && result <= 1.0);
    return result;
}

double Candidates::compute_noncandidate_conservative_stepsize(
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> displacements,
    const double dhat) const
{
    assert(displacements.rows() == mesh.num_vertices());

    if (empty()) {
        return 1; // No possible collisions, so can take full step.
    }

    const Eigen::MatrixXi& E = mesh.edges();
    const Eigen::MatrixXi& F = mesh.faces();

    std::vector<std::atomic<bool>> is_vertex_a_candidates(mesh.num_vertices());
    for (size_t i = 0; i < mesh.num_vertices(); ++i) {
        is_vertex_a_candidates[i].store(false, std::memory_order_relaxed);
    }

    tbb::parallel_for(size_t(0), size(), [&](size_t i) {
        for (const index_t vid : (*this)[i].vertex_ids(E, F)) {
            if (vid < 0) {
                break;
            }
            is_vertex_a_candidates[vid].store(true, std::memory_order_relaxed);
        }
    });

    double max_displacement = tbb::parallel_reduce(
        tbb::blocked_range<size_t>(0, displacements.rows()), 0.0,
        [&](const tbb::blocked_range<size_t>& r, double local_max) -> double {
            for (size_t i = r.begin(); i < r.end(); ++i) {
                if (!is_vertex_a_candidates[i]) {
                    continue;
                }
                local_max = std::max(local_max, displacements.row(i).norm());
            }
            return local_max;
        },
        [](double a, double b) { return std::max(a, b); });

    return 0.5 * dhat / max_displacement;
}

double Candidates::compute_cfl_stepsize(
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices_t0,
    Eigen::ConstRef<Eigen::MatrixXd> vertices_t1,
    const double dhat,
    const double min_distance,
    BroadPhase* broad_phase,
    const NarrowPhaseCCD& narrow_phase_ccd) const
{
    assert(vertices_t0.rows() == mesh.num_vertices());
    assert(vertices_t1.rows() == mesh.num_vertices());

    const double alpha_C = this->compute_collision_free_stepsize(
        mesh, vertices_t0, vertices_t1, min_distance, narrow_phase_ccd);

    const double alpha_F = this->compute_noncandidate_conservative_stepsize(
        mesh, vertices_t1 - vertices_t0, dhat);

    // If alpha_F < 0.5 * alpha_C, then we should do full CCD.
    if (alpha_F < 0.5 * alpha_C) {
        return ipc::compute_collision_free_stepsize(
            mesh, vertices_t0, vertices_t1, min_distance, broad_phase,
            narrow_phase_ccd);
    }
    return std::min(alpha_C, alpha_F);
}

Eigen::VectorXd Candidates::compute_per_vertex_safe_distances(
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices,
    const double inflation_radius,
    const double min_distance) const
{
    assert(inflation_radius >= min_distance);

    // Initialize atomic distances for each vertex
    std::vector<std::atomic<double>> min_distances(mesh.num_vertices());
    for (size_t i = 0; i < mesh.num_vertices(); ++i) {
        min_distances[i].store(
            inflation_radius - min_distance, std::memory_order_relaxed);
    }

    tbb::parallel_for(size_t(0), size(), [&](size_t i) {
        const CollisionStencil& candidate = (*this)[i];

        const double d = sqrt(candidate.compute_distance(
                             vertices, mesh.edges(), mesh.faces()))
            - min_distance;

        // Compute the distance for each vertex in the candidate
        for (auto vid : candidate.vertex_ids(mesh.edges(), mesh.faces())) {
            if (vid < 0) {
                break; // No more vertices in this candidate
            }
            // Update the minimum distance atomically
            double old_val = min_distances[vid].load(std::memory_order_relaxed);
            while (d < old_val
                   && !min_distances[vid].compare_exchange_weak(
                       old_val, d, std::memory_order_relaxed)) { }
        }
    });

    // Convert atomic distances to a vector
    Eigen::VectorXd result(mesh.num_vertices());
    for (size_t i = 0; i < mesh.num_vertices(); ++i) {
        result[i] = min_distances[i].load(std::memory_order_relaxed);
    }
    assert((result.array() >= 0).all());

    return result;
}

// Eigen::VectorXd Candidates::compute_per_vertex_collision_free_stepsize(
//     const CollisionMesh& mesh,
//     Eigen::ConstRef<Eigen::MatrixXd> vertices_t0,
//     Eigen::ConstRef<Eigen::MatrixXd> vertices_t1,
//     const double min_distance,
//     const NarrowPhaseCCD& narrow_phase_ccd) const
// {
//     // Initialize atomic step size for each vertex
//     std::vector<std::atomic<double>> max_step(mesh.num_vertices());
//     for (size_t i = 0; i < mesh.num_vertices(); ++i) {
//         max_step[i].store(1.0, std::memory_order_relaxed);
//     }

//     tbb::parallel_for(
//         tbb::blocked_range<size_t>(0, size()),
//         [&](const tbb::blocked_range<size_t>& r) {
//             for (size_t i = r.begin(); i < r.end(); i++) {
//                 const CollisionStencil& candidate = (*this)[i];
//                 const auto vertex_ids =
//                     candidate.vertex_ids(mesh.edges(), mesh.faces());

//                 double tmax = 1.0;
//                 for (const index_t vid : vertex_ids) {
//                     if (vid < 0) {
//                         break; // No more vertices in this candidate
//                     }
//                     // Get the maximum time of impact for this vertex
//                     tmax = std::min(
//                         tmax, max_step[vid].load(std::memory_order_relaxed));
//                 }

//                 double toi;
//                 const bool collides = candidate.ccd(
//                     candidate.dof(vertices_t0, mesh.edges(), mesh.faces()),
//                     candidate.dof(vertices_t1, mesh.edges(), mesh.faces()),
//                     toi, min_distance, tmax, narrow_phase_ccd);

//                 if (collides) {
//                     // Compute the distance for each vertex in the candidate
//                     for (const index_t vid : vertex_ids) {
//                         if (vid < 0) {
//                             break; // No more vertices in this candidate
//                         }
//                         // Update the max_step atomically
//                         double old_val =
//                             max_step[vid].load(std::memory_order_relaxed);
//                         while (toi < old_val
//                                && !max_step[vid].compare_exchange_weak(
//                                    old_val, toi, std::memory_order_relaxed))
//                                    { }
//                     }
//                 }
//             }
//         });

//     // Convert atomic distances to a vector
//     Eigen::VectorXd result(mesh.num_vertices());
//     for (size_t i = 0; i < mesh.num_vertices(); ++i) {
//         result[i] = max_step[i].load(std::memory_order_relaxed);
//     }
//     assert((result.array() >= 0).all());

//     return result;
// }

// ============================================================================

size_t Candidates::size() const
{
    return vv_candidates.size() + ev_candidates.size() + ee_candidates.size()
        + fv_candidates.size() + pv_candidates.size();
}

bool Candidates::empty() const
{
    return vv_candidates.empty() && ev_candidates.empty()
        && ee_candidates.empty() && fv_candidates.empty()
        && pv_candidates.empty();
}

void Candidates::clear()
{
    vv_candidates.clear();
    ev_candidates.clear();
    ee_candidates.clear();
    fv_candidates.clear();
    ef_candidates.clear();
    ff_candidates.clear();
    pv_candidates.clear();
}

CollisionStencil& Candidates::operator[](size_t i)
{
    if (i < vv_candidates.size()) {
        return vv_candidates[i];
    }
    i -= vv_candidates.size();
    if (i < ev_candidates.size()) {
        return ev_candidates[i];
    }
    i -= ev_candidates.size();
    if (i < ee_candidates.size()) {
        return ee_candidates[i];
    }
    i -= ee_candidates.size();
    if (i < fv_candidates.size()) {
        return fv_candidates[i];
    }
    i -= fv_candidates.size();
    if (i < pv_candidates.size()) {
        return pv_candidates[i];
    }
    throw std::out_of_range("Candidate index is out of range!");
}

const CollisionStencil& Candidates::operator[](size_t i) const
{
    if (i < vv_candidates.size()) {
        return vv_candidates[i];
    }
    i -= vv_candidates.size();
    if (i < ev_candidates.size()) {
        return ev_candidates[i];
    }
    i -= ev_candidates.size();
    if (i < ee_candidates.size()) {
        return ee_candidates[i];
    }
    i -= ee_candidates.size();
    if (i < fv_candidates.size()) {
        return fv_candidates[i];
    }
    i -= fv_candidates.size();
    if (i < pv_candidates.size()) {
        return pv_candidates[i];
    }
    throw std::out_of_range("Candidate index is out of range!");
}

bool Candidates::is_vertex_vertex(size_t i) const
{
    return i < vv_candidates.size();
}

bool Candidates::is_edge_vertex(size_t i) const
{
    return i >= vv_candidates.size()
        && i < vv_candidates.size() + ev_candidates.size();
}

bool Candidates::is_edge_edge(size_t i) const
{
    return i >= vv_candidates.size() + ev_candidates.size()
        && i
        < vv_candidates.size() + ev_candidates.size() + ee_candidates.size();
}

bool Candidates::is_face_vertex(size_t i) const
{
    return i
        >= vv_candidates.size() + ev_candidates.size() + ee_candidates.size()
        && i < vv_candidates.size() + ev_candidates.size()
            + ee_candidates.size() + fv_candidates.size();
}

bool Candidates::is_plane_vertex(size_t i) const
{
    return i >= vv_candidates.size() + ev_candidates.size()
            + ee_candidates.size() + fv_candidates.size()
        && i < vv_candidates.size() + ev_candidates.size()
            + ee_candidates.size() + fv_candidates.size()
            + pv_candidates.size();
}

// == Convert to subelement candidates ========================================

namespace {

    /// @brief Convert element-vertex candidates to vertex-vertex candidates
    /// @param elements Elements matrix of the mesh
    /// @param vertices Vertex positions of the mesh
    /// @param ev_candidates Element-vertex candidates
    /// @param is_active Function to determine if a candidate is active
    /// @return Vertex-vertex candidates
    template <typename Candidate>
    std::vector<VertexVertexCandidate>
    element_vertex_to_vertex_vertex_candidates(
        Eigen::ConstRef<Eigen::MatrixXi> elements,
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        const std::vector<Candidate>& candidates,
        const std::function<bool(double)>& is_active)
    {
        std::vector<VertexVertexCandidate> vv_candidates;
        for (const auto& [ei, vi] : candidates) {
            for (int j = 0; j < elements.cols(); j++) {
                const int vj = elements(ei, j);
                if (is_active(point_point_distance(
                        vertices.row(vi), vertices.row(vj)))) {
                    vv_candidates.emplace_back(vi, vj);
                }
            }
        }

        // Remove duplicates
        tbb::parallel_sort(vv_candidates.begin(), vv_candidates.end());
        vv_candidates.erase(
            std::unique(vv_candidates.begin(), vv_candidates.end()),
            vv_candidates.end());

        return vv_candidates;
    }

} // namespace

std::vector<VertexVertexCandidate> Candidates::edge_vertex_to_vertex_vertex(
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices,
    const std::function<bool(double)>& is_active) const
{
    return element_vertex_to_vertex_vertex_candidates(
        mesh.edges(), vertices, ev_candidates, is_active);
}

std::vector<VertexVertexCandidate> Candidates::face_vertex_to_vertex_vertex(
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices,
    const std::function<bool(double)>& is_active) const
{
    return element_vertex_to_vertex_vertex_candidates(
        mesh.faces(), vertices, fv_candidates, is_active);
}

std::vector<EdgeVertexCandidate> Candidates::face_vertex_to_edge_vertex(
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices,
    const std::function<bool(double)>& is_active) const
{
    std::vector<EdgeVertexCandidate> C_ev;
    for (const auto& [fi, vi] : fv_candidates) {
        for (int j = 0; j < 3; j++) {
            const int ei = mesh.faces_to_edges()(fi, j);
            const int vj = mesh.edges()(ei, 0);
            const int vk = mesh.edges()(ei, 1);
            if (is_active(point_edge_distance(
                    vertices.row(vi), vertices.row(vj), vertices.row(vk)))) {
                C_ev.emplace_back(ei, vi);
            }
        }
    }

    // Remove duplicates
    tbb::parallel_sort(C_ev.begin(), C_ev.end());
    C_ev.erase(std::unique(C_ev.begin(), C_ev.end()), C_ev.end());

    return C_ev;
}

std::vector<EdgeVertexCandidate> Candidates::edge_edge_to_edge_vertex(
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices,
    const std::function<bool(double)>& is_active) const
{
    std::vector<EdgeVertexCandidate> C_ev;
    for (const EdgeEdgeCandidate& ee : ee_candidates) {
        for (int i = 0; i < 2; i++) {
            const int ei = i == 0 ? ee.edge0_id : ee.edge1_id;
            const int ej = i == 0 ? ee.edge1_id : ee.edge0_id;

            const int ei0 = mesh.edges()(ei, 0);
            const int ei1 = mesh.edges()(ei, 1);

            for (int j = 0; j < 2; j++) {
                const int vj = mesh.edges()(ej, j);
                if (is_active(point_edge_distance(
                        vertices.row(vj), vertices.row(ei0),
                        vertices.row(ei1)))) {
                    C_ev.emplace_back(ei, vj);
                }
            }
        }
    }

    // Remove duplicates
    tbb::parallel_sort(C_ev.begin(), C_ev.end());
    C_ev.erase(std::unique(C_ev.begin(), C_ev.end()), C_ev.end());

    return C_ev;
}

// ============================================================================

bool Candidates::save_obj(
    const std::string& filename,
    Eigen::ConstRef<Eigen::MatrixXd> vertices,
    Eigen::ConstRef<Eigen::MatrixXi> edges,
    Eigen::ConstRef<Eigen::MatrixXi> faces) const
{
    std::ofstream obj(filename, std::ios::out);
    if (!obj.is_open()) {
        return false;
    }
    int v_offset = 0;
    ipc::save_obj(obj, vertices, edges, faces, vv_candidates, v_offset);
    v_offset += vv_candidates.size() * 2;
    ipc::save_obj(obj, vertices, edges, faces, ev_candidates, v_offset);
    v_offset += ev_candidates.size() * 3;
    ipc::save_obj(obj, vertices, edges, faces, ee_candidates, v_offset);
    v_offset += ee_candidates.size() * 4;
    ipc::save_obj(obj, vertices, faces, faces, fv_candidates, v_offset);
    return true;
}

void Candidates::convert_candidates_to_sets()
{
    m_sets = std::make_shared<AdjacencySets>();

    for (const auto& vv : vv_candidates) {
        m_sets->vv[vv.vertex0_id].insert(vv.vertex1_id);
        m_sets->vv[vv.vertex1_id].insert(vv.vertex0_id);
    }
    for (const auto& ee : ee_candidates) {
        m_sets->ee[ee.edge0_id].insert(ee.edge1_id);
        m_sets->ee[ee.edge1_id].insert(ee.edge0_id);
    }
    for (const auto& ff : ff_candidates) {
        m_sets->ff[ff.face0_id].insert(ff.face1_id);
        m_sets->ff[ff.face1_id].insert(ff.face0_id);
    }
    for (const auto& ev : ev_candidates) {
        m_sets->ev[ev.edge_id].insert(ev.vertex_id);
        m_sets->ve[ev.vertex_id].insert(ev.edge_id);
    }
    for (const auto& fv : fv_candidates) {
        m_sets->fv[fv.face_id].insert(fv.vertex_id);
        m_sets->vf[fv.vertex_id].insert(fv.face_id);
    }
    for (const auto& ef : ef_candidates) {
        m_sets->ef[ef.edge_id].insert(ef.face_id);
        m_sets->fe[ef.face_id].insert(ef.edge_id);
    }
}

std::set<index_t> Candidates::vv_set(index_t id) const
{
    if (!m_sets) {
        return {};
    }

    assert(m_mesh.num_vertices());
    std::set<index_t> out;
    if (auto iter = m_sets->vv.find(id); iter != m_sets->vv.end()) {
        out = iter->second;
    }

    if (m_mesh.dim() == 2) {
        for (const index_t ej : ve_set(id)) {
            out.insert(m_mesh.edges()(ej, 0));
            out.insert(m_mesh.edges()(ej, 1));
        }
    }
    out.erase(id);
    return out;
}
std::set<index_t> Candidates::ve_set(index_t id) const
{
    if (!m_sets) {
        return {};
    }

    if (auto iter = m_sets->ve.find(id); iter != m_sets->ve.end()) {
        return iter->second;
    }
    return {};
}
std::set<index_t> Candidates::vf_set(index_t id) const
{
    if (!m_sets) {
        return {};
    }

    if (auto iter = m_sets->vf.find(id); iter != m_sets->vf.end()) {
        return iter->second;
    }
    return {};
}

std::set<index_t> Candidates::ev_set(index_t id) const
{
    if (!m_sets) {
        return {};
    }

    assert(m_mesh.num_vertices());
    std::set<index_t> out;
    if (auto iter = m_sets->ev.find(id); iter != m_sets->ev.end()) {
        out = iter->second;
    }
    for (index_t lv = 0; lv < 2; ++lv) {
        out.insert(m_mesh.edges()(id, lv));
    }
    return out;
}
std::set<index_t> Candidates::ee_set(index_t id) const
{
    if (!m_sets) {
        return {};
    }

    assert(m_mesh.num_vertices());
    std::set<index_t> out;
    if (auto iter = m_sets->ee.find(id); iter != m_sets->ee.end()) {
        out = iter->second;
    }
    for (index_t lv = 0; lv < 2; ++lv) {
        for (index_t eid : m_mesh.vertices_to_edges()[m_mesh.edges()(id, lv)]) {
            out.insert(eid);
        }
    }
    // In 2D, EE candidates are never built by the broad phase. Reconstruct
    // them from EV candidates symmetrically:
    // (a) edges adjacent to vertices that are close to edge id (via ev_set)
    // (b) edges that id's own endpoints are close to (via ve_set)
    if (m_mesh.dim() == 2) {
        for (const index_t vj : ev_set(id)) {
            for (const index_t ej : m_mesh.vertices_to_edges()[vj]) {
                out.insert(ej);
            }
        }
        for (index_t lv = 0; lv < 2; ++lv) {
            const index_t vi = m_mesh.edges()(id, lv);
            for (const index_t ej : ve_set(vi)) {
                out.insert(ej);
            }
        }
    }
    out.erase(id);
    return out;
}
std::set<index_t> Candidates::ef_set(index_t id) const
{
    if (!m_sets) {
        return {};
    }

    assert(m_mesh.num_vertices());
    std::set<index_t> out;
    if (auto iter = m_sets->ef.find(id); iter != m_sets->ef.end()) {
        out = iter->second;
    }
    for (index_t lv = 0; lv < 2; ++lv) {
        const auto& faces = m_mesh.vertices_to_faces()[m_mesh.edges()(id, lv)];
        for (int fid : faces) {
            out.insert(fid);
        }
    }
    for (const index_t fid : m_mesh.edges_to_faces()[id]) {
        out.erase(fid);
    }
    return out;
}

std::set<index_t> Candidates::fv_set(index_t id) const
{
    if (!m_sets) {
        return {};
    }

    assert(m_mesh.num_vertices());
    std::set<index_t> out;
    if (auto iter = m_sets->fv.find(id); iter != m_sets->fv.end()) {
        out = iter->second;
    }
    for (index_t lv = 0; lv < 3; ++lv) {
        out.insert(m_mesh.faces()(id, lv));
    }
    return out;
}
std::set<index_t> Candidates::fe_set(index_t id) const
{
    if (!m_sets) {
        return {};
    }

    assert(m_mesh.num_vertices());
    std::set<index_t> out;
    if (auto iter = m_sets->fe.find(id); iter != m_sets->fe.end()) {
        out = iter->second;
    }
    for (index_t lv = 0; lv < 3; ++lv) {
        for (index_t eid : m_mesh.vertices_to_edges()[m_mesh.faces()(id, lv)]) {
            out.insert(eid);
        }
    }
    return out;
}
std::set<index_t> Candidates::ff_set(index_t id) const
{
    if (!m_sets) {
        return {};
    }

    assert(m_mesh.num_vertices());
    std::set<index_t> out;
    if (auto iter = m_sets->ff.find(id); iter != m_sets->ff.end()) {
        out = iter->second;
    }
    for (index_t lv = 0; lv < 3; ++lv) {
        const index_t vid = m_mesh.faces()(id, lv);
        for (index_t fid : m_mesh.vertices_to_faces()[vid]) {
            out.insert(fid);
        }
    }
    out.erase(id);
    return out;
}

} // namespace ipc
