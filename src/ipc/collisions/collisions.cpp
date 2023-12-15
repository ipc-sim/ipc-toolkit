#include "collisions.hpp"

#include <ipc/collisions/collisions_builder.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/distance/point_line.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/point_plane.hpp>
#include <ipc/utils/local_to_global.hpp>

#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>
#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>

#include <stdexcept> // std::out_of_range

namespace ipc {

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
        const Eigen::MatrixXi& elements,
        const Eigen::MatrixXd& vertices,
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

    std::vector<VertexVertexCandidate> edge_vertex_to_vertex_vertex_candidates(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const std::vector<EdgeVertexCandidate>& ev_candidates,
        const std::function<bool(double)>& is_active)
    {
        return element_vertex_to_vertex_vertex_candidates(
            mesh.edges(), vertices, ev_candidates, is_active);
    }

    std::vector<VertexVertexCandidate> face_vertex_to_vertex_vertex_candidates(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const std::vector<FaceVertexCandidate>& fv_candidates,
        const std::function<bool(double)>& is_active)
    {
        return element_vertex_to_vertex_vertex_candidates(
            mesh.faces(), vertices, fv_candidates, is_active);
    }

    std::vector<EdgeVertexCandidate> face_vertex_to_edge_vertex_candidates(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const std::vector<FaceVertexCandidate>& fv_candidates,
        const std::function<bool(double)>& is_active)
    {
        std::vector<EdgeVertexCandidate> ev_candidates;
        for (const auto& [fi, vi] : fv_candidates) {
            for (int j = 0; j < 3; j++) {
                const int ei = mesh.faces_to_edges()(fi, j);
                const int vj = mesh.edges()(ei, 0);
                const int vk = mesh.edges()(ei, 1);
                if (is_active(point_edge_distance(
                        vertices.row(vi), //
                        vertices.row(vj), vertices.row(vk)))) {
                    ev_candidates.emplace_back(ei, vi);
                }
            }
        }

        // Remove duplicates
        tbb::parallel_sort(ev_candidates.begin(), ev_candidates.end());
        ev_candidates.erase(
            std::unique(ev_candidates.begin(), ev_candidates.end()),
            ev_candidates.end());

        return ev_candidates;
    }

    std::vector<EdgeVertexCandidate> edge_edge_to_edge_vertex_candidates(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const std::vector<EdgeEdgeCandidate>& ee_candidates,
        const std::function<bool(double)>& is_active)
    {
        std::vector<EdgeVertexCandidate> ev_candidates;
        for (const EdgeEdgeCandidate& ee : ee_candidates) {
            for (int i = 0; i < 2; i++) {
                const int ei = i == 0 ? ee.edge0_id : ee.edge1_id;
                const int ej = i == 0 ? ee.edge1_id : ee.edge0_id;

                const int ei0 = mesh.edges()(ei, 0);
                const int ei1 = mesh.edges()(ei, 1);

                for (int j = 0; j < 2; j++) {
                    const int vj = mesh.edges()(ej, j);
                    if (is_active(point_edge_distance(
                            vertices.row(vj), //
                            vertices.row(ei0), vertices.row(ei1)))) {
                        ev_candidates.emplace_back(ei, vj);
                    }
                }
            }
        }

        // Remove duplicates
        tbb::parallel_sort(ev_candidates.begin(), ev_candidates.end());
        ev_candidates.erase(
            std::unique(ev_candidates.begin(), ev_candidates.end()),
            ev_candidates.end());

        return ev_candidates;
    }
} // namespace

void Collisions::build(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const double dhat,
    const double dmin,
    const BroadPhaseMethod broad_phase_method)
{
    assert(vertices.rows() == mesh.num_vertices());

    double inflation_radius = (dhat + dmin) / 2;

    Candidates candidates;
    candidates.build(mesh, vertices, inflation_radius, broad_phase_method);

    this->build(candidates, mesh, vertices, dhat, dmin);
}

void Collisions::build(
    const Candidates& candidates,
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const double dhat,
    const double dmin)
{
    assert(vertices.rows() == mesh.num_vertices());

    clear();

    // Cull the candidates by measuring the distance and dropping those that are
    // greater than dhat.
    const double offset_sqr = (dmin + dhat) * (dmin + dhat);
    auto is_active = [&](double distance_sqr) {
        return distance_sqr < offset_sqr;
    };

    tbb::enumerable_thread_specific<CollisionsBuilder> storage(
        use_convergent_formulation(), are_shape_derivatives_enabled());

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), candidates.vv_candidates.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            storage.local().add_vertex_vertex_collisions(
                mesh, vertices, candidates.vv_candidates, is_active, r.begin(),
                r.end());
        });

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), candidates.ev_candidates.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            storage.local().add_edge_vertex_collisions(
                mesh, vertices, candidates.ev_candidates, is_active, r.begin(),
                r.end());
        });

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), candidates.ee_candidates.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            storage.local().add_edge_edge_collisions(
                mesh, vertices, candidates.ee_candidates, is_active, r.begin(),
                r.end());
        });

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), candidates.fv_candidates.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            storage.local().add_face_vertex_collisions(
                mesh, vertices, candidates.fv_candidates, is_active, r.begin(),
                r.end());
        });

    if (use_convergent_formulation()) {
        if (candidates.ev_candidates.size() > 0) {
            // Convert edge-vertex to vertex-vertex
            const std::vector<VertexVertexCandidate> vv_candidates =
                edge_vertex_to_vertex_vertex_candidates(
                    mesh, vertices, candidates.ev_candidates, is_active);

            tbb::parallel_for(
                tbb::blocked_range<size_t>(size_t(0), vv_candidates.size()),
                [&](const tbb::blocked_range<size_t>& r) {
                    storage.local()
                        .add_edge_vertex_negative_vertex_vertex_collisions(
                            mesh, vertices, vv_candidates, r.begin(), r.end());
                });
        }

        if (candidates.ee_candidates.size() > 0) {
            // Convert edge-edge to edge-vertex
            const auto ev_candidates = edge_edge_to_edge_vertex_candidates(
                mesh, vertices, candidates.ee_candidates, is_active);

            tbb::parallel_for(
                tbb::blocked_range<size_t>(size_t(0), ev_candidates.size()),
                [&](const tbb::blocked_range<size_t>& r) {
                    storage.local()
                        .add_edge_edge_negative_edge_vertex_collisions(
                            mesh, vertices, ev_candidates, r.begin(), r.end());
                });
        }

        if (candidates.fv_candidates.size() > 0) {
            // Convert face-vertex to edge-vertex
            const std::vector<EdgeVertexCandidate> ev_candidates =
                face_vertex_to_edge_vertex_candidates(
                    mesh, vertices, candidates.fv_candidates, is_active);

            tbb::parallel_for(
                tbb::blocked_range<size_t>(size_t(0), ev_candidates.size()),
                [&](const tbb::blocked_range<size_t>& r) {
                    storage.local()
                        .add_face_vertex_negative_edge_vertex_collisions(
                            mesh, vertices, ev_candidates, r.begin(), r.end());
                });

            // Convert face-vertex to vertex-vertex
            const std::vector<VertexVertexCandidate> vv_candidates =
                face_vertex_to_vertex_vertex_candidates(
                    mesh, vertices, candidates.fv_candidates, is_active);

            tbb::parallel_for(
                tbb::blocked_range<size_t>(size_t(0), vv_candidates.size()),
                [&](const tbb::blocked_range<size_t>& r) {
                    storage.local()
                        .add_face_vertex_positive_vertex_vertex_collisions(
                            mesh, vertices, vv_candidates, r.begin(), r.end());
                });
        }
    }

    // -------------------------------------------------------------------------

    CollisionsBuilder::merge(storage, *this);

    // logger().debug(to_string(mesh, vertices));

    for (size_t ci = 0; ci < size(); ci++) {
        Collision& collision = (*this)[ci];
        collision.dmin = dmin;
    }

    if (use_convergent_formulation()) {
        // NOTE: When using the convergent formulation we want the barrier to
        // have units of Pa⋅m, so κ gets units of Pa and the barrier function
        // should have units of m. See notebooks/physical_barrier.ipynb for more
        // details.
        const double barrier_to_physical_barrier_divisor =
            dhat * std::pow(dhat + 2 * dmin, 2);

        for (size_t ci = 0; ci < size(); ci++) {
            Collision& collision = (*this)[ci];
            collision.weight /= barrier_to_physical_barrier_divisor;
            if (are_shape_derivatives_enabled()) {
                collision.weight_gradient /=
                    barrier_to_physical_barrier_divisor;
            }
        }
    }
}

void Collisions::set_use_convergent_formulation(
    const bool use_convergent_formulation)
{
    if (!empty()
        && use_convergent_formulation != m_use_convergent_formulation) {
        logger().warn(
            "Setting use_convergent_formulation after building collisions. "
            "Re-build collisions for this to have an effect.");
    }

    m_use_convergent_formulation = use_convergent_formulation;
}

void Collisions::set_are_shape_derivatives_enabled(
    const bool are_shape_derivatives_enabled)
{
    if (!empty()
        && are_shape_derivatives_enabled != m_are_shape_derivatives_enabled) {
        logger().warn(
            "Setting enable_shape_derivatives after building collisions. "
            "Re-build collisions for this to have an effect.");
    }

    m_are_shape_derivatives_enabled = are_shape_derivatives_enabled;
}

// ============================================================================

// NOTE: Actually distance squared
double Collisions::compute_minimum_distance(
    const CollisionMesh& mesh, const Eigen::MatrixXd& vertices) const
{
    assert(vertices.rows() == mesh.num_vertices());

    if (empty()) {
        return std::numeric_limits<double>::infinity();
    }

    const Eigen::MatrixXi& edges = mesh.edges();
    const Eigen::MatrixXi& faces = mesh.faces();

    tbb::enumerable_thread_specific<double> storage(
        std::numeric_limits<double>::infinity());

    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, size()),
        [&](tbb::blocked_range<size_t> r) {
            double& local_min_dist = storage.local();

            for (size_t i = r.begin(); i < r.end(); i++) {
                const double dist = (*this)[i].compute_distance(
                    (*this)[i].dof(vertices, edges, faces));

                if (dist < local_min_dist) {
                    local_min_dist = dist;
                }
            }
        });

    return storage.combine([](double a, double b) { return std::min(a, b); });
}

// ============================================================================

size_t Collisions::size() const
{
    return vv_collisions.size() + ev_collisions.size() + ee_collisions.size()
        + fv_collisions.size() + pv_collisions.size();
}

bool Collisions::empty() const
{
    return vv_collisions.empty() && ev_collisions.empty()
        && ee_collisions.empty() && fv_collisions.empty()
        && pv_collisions.empty();
}

void Collisions::clear()
{
    vv_collisions.clear();
    ev_collisions.clear();
    ee_collisions.clear();
    fv_collisions.clear();
    pv_collisions.clear();
}

Collision& Collisions::operator[](size_t i)
{
    if (i < vv_collisions.size()) {
        return vv_collisions[i];
    }
    i -= vv_collisions.size();
    if (i < ev_collisions.size()) {
        return ev_collisions[i];
    }
    i -= ev_collisions.size();
    if (i < ee_collisions.size()) {
        return ee_collisions[i];
    }
    i -= ee_collisions.size();
    if (i < fv_collisions.size()) {
        return fv_collisions[i];
    }
    i -= fv_collisions.size();
    if (i < pv_collisions.size()) {
        return pv_collisions[i];
    }
    throw std::out_of_range("Collision index is out of range!");
}

const Collision& Collisions::operator[](size_t i) const
{
    if (i < vv_collisions.size()) {
        return vv_collisions[i];
    }
    i -= vv_collisions.size();
    if (i < ev_collisions.size()) {
        return ev_collisions[i];
    }
    i -= ev_collisions.size();
    if (i < ee_collisions.size()) {
        return ee_collisions[i];
    }
    i -= ee_collisions.size();
    if (i < fv_collisions.size()) {
        return fv_collisions[i];
    }
    i -= fv_collisions.size();
    if (i < pv_collisions.size()) {
        return pv_collisions[i];
    }
    throw std::out_of_range("Collision index is out of range!");
}

bool Collisions::is_vertex_vertex(size_t i) const
{
    return i < vv_collisions.size();
}

bool Collisions::is_edge_vertex(size_t i) const
{
    return i >= vv_collisions.size()
        && i < vv_collisions.size() + ev_collisions.size();
}

bool Collisions::is_edge_edge(size_t i) const
{
    return i >= vv_collisions.size() + ev_collisions.size()
        && i
        < vv_collisions.size() + ev_collisions.size() + ee_collisions.size();
}

bool Collisions::is_face_vertex(size_t i) const
{
    return i
        >= vv_collisions.size() + ev_collisions.size() + ee_collisions.size()
        && i < vv_collisions.size() + ev_collisions.size()
            + ee_collisions.size() + fv_collisions.size();
}

bool Collisions::is_plane_vertex(size_t i) const
{
    return i >= vv_collisions.size() + ev_collisions.size()
            + ee_collisions.size() + fv_collisions.size()
        && i < vv_collisions.size() + ev_collisions.size()
            + ee_collisions.size() + fv_collisions.size()
            + pv_collisions.size();
}

std::string Collisions::to_string(
    const CollisionMesh& mesh, const Eigen::MatrixXd& vertices) const
{
    std::stringstream ss;
    for (const auto& vv : vv_collisions) {
        ss << "\n"
           << fmt::format(
                  "vv: {} {}, w: {:g}, d: {:g}", vv.vertex0_id, vv.vertex1_id,
                  vv.weight,
                  vv.compute_distance(
                      vv.dof(vertices, mesh.edges(), mesh.faces())));
    }
    for (const auto& ev : ev_collisions) {
        ss << "\n"
           << fmt::format(
                  "ev: {}=({}, {}) {}, w: {:g}, d: {:g}", ev.edge_id,
                  mesh.edges()(ev.edge_id, 0), mesh.edges()(ev.edge_id, 1),
                  ev.vertex_id, ev.weight,
                  ev.compute_distance(
                      ev.dof(vertices, mesh.edges(), mesh.faces())));
    }
    for (const auto& ee : ee_collisions) {
        ss << "\n"
           << fmt::format(
                  "ee: {}=({}, {}) {}=({}, {}), w: {:g}, dtype: {}, d: {:g}",
                  ee.edge0_id, mesh.edges()(ee.edge0_id, 0),
                  mesh.edges()(ee.edge0_id, 1), ee.edge1_id,
                  mesh.edges()(ee.edge1_id, 0), mesh.edges()(ee.edge1_id, 1),
                  ee.weight, int(ee.dtype),
                  ee.compute_distance(
                      ee.dof(vertices, mesh.edges(), mesh.faces())));
    }
    for (const auto& fv : fv_collisions) {
        ss << "\n"
           << fmt::format(
                  "fv: {}=({}, {}, {}) {}, w: {:g}, d: {:g}", fv.face_id,
                  mesh.faces()(fv.face_id, 0), mesh.faces()(fv.face_id, 1),
                  mesh.faces()(fv.face_id, 2), fv.vertex_id, fv.weight,
                  fv.compute_distance(
                      fv.dof(vertices, mesh.edges(), mesh.faces())));
    }
    return ss.str();
}

} // namespace ipc
