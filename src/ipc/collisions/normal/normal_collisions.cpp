#include "normal_collisions.hpp"

#include <ipc/collisions/normal/normal_collisions_builder.hpp>
#include <ipc/utils/local_to_global.hpp>

#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>

#include <stdexcept> // std::out_of_range

namespace ipc {

namespace {
    inline double sqr(double x) { return x * x; }
} // namespace

void NormalCollisions::build(
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices,
    const double dhat,
    const double dmin,
    BroadPhase* broad_phase)
{
    assert(vertices.rows() == mesh.num_vertices());

    const double inflation_radius = 0.5 * (dhat + dmin);

    Candidates candidates;
    candidates.build(mesh, vertices, inflation_radius, broad_phase);

    this->build(candidates, mesh, vertices, dhat, dmin);
}

void NormalCollisions::build(
    const Candidates& candidates,
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices,
    const double dhat,
    const double dmin)
{
    assert(vertices.rows() == mesh.num_vertices());

    clear();

    // Cull the candidates by measuring the distance and dropping those that are
    // greater than dhat.
    auto is_active = [offset_sqr = sqr(dmin + dhat)](double distance_sqr) {
        return distance_sqr < offset_sqr;
    };

    tbb::enumerable_thread_specific<NormalCollisionsBuilder> storage(
        use_area_weighting(), enable_shape_derivatives(),
        collision_set_type() == CollisionSetType::OGC);

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

    if (collision_set_type() == CollisionSetType::IMPROVED_MAX_APPROX) {
        if (!candidates.ev_candidates.empty()) {
            // Convert edge-vertex to vertex-vertex
            const auto vv_candidates = candidates.edge_vertex_to_vertex_vertex(
                mesh, vertices, is_active);

            tbb::parallel_for(
                tbb::blocked_range<size_t>(size_t(0), vv_candidates.size()),
                [&](const tbb::blocked_range<size_t>& r) {
                    storage.local()
                        .add_edge_vertex_negative_vertex_vertex_collisions(
                            mesh, vertices, vv_candidates, r.begin(), r.end());
                });
        }

        if (!candidates.ee_candidates.empty()) {
            // Convert edge-edge to edge-vertex
            const auto ev_candidates =
                candidates.edge_edge_to_edge_vertex(mesh, vertices, is_active);

            tbb::parallel_for(
                tbb::blocked_range<size_t>(size_t(0), ev_candidates.size()),
                [&](const tbb::blocked_range<size_t>& r) {
                    storage.local()
                        .add_edge_edge_negative_edge_vertex_collisions(
                            mesh, vertices, ev_candidates, r.begin(), r.end());
                });
        }

        if (!candidates.fv_candidates.empty()) {
            // Convert face-vertex to edge-vertex
            const auto ev_candidates = candidates.face_vertex_to_edge_vertex(
                mesh, vertices, is_active);

            tbb::parallel_for(
                tbb::blocked_range<size_t>(size_t(0), ev_candidates.size()),
                [&](const tbb::blocked_range<size_t>& r) {
                    storage.local()
                        .add_face_vertex_negative_edge_vertex_collisions(
                            mesh, vertices, ev_candidates, r.begin(), r.end());
                });

            // Convert face-vertex to vertex-vertex
            const auto vv_candidates = candidates.face_vertex_to_vertex_vertex(
                mesh, vertices, is_active);

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

    NormalCollisionsBuilder::merge(storage, *this);

    // logger().debug(to_string(mesh, vertices));

    for (size_t ci = 0; ci < size(); ci++) {
        NormalCollision& collision = (*this)[ci];
        collision.dmin = dmin;
    }
}

void NormalCollisions::set_use_area_weighting(const bool use_area_weighting)
{
    if (!empty() && use_area_weighting != m_use_area_weighting) {
        logger().warn(
            "Setting use_area_weighting after building collisions. "
            "Re-build collisions for this to have an effect.");
    }

    if (!use_area_weighting
        && collision_set_type() == CollisionSetType::IMPROVED_MAX_APPROX) {
        logger().warn(
            "Disabling area weighting while using the improved max approximation may lead to incorrect results.");
    }

    m_use_area_weighting = use_area_weighting;
}

void NormalCollisions::set_collision_set_type(const CollisionSetType type)
{
    if (!empty() && type != m_collision_set_type) {
        logger().warn(
            "Setting collision_set_type after building collisions. "
            "Re-build collisions for this to have an effect.");
    }

    if (!use_area_weighting()
        && type == CollisionSetType::IMPROVED_MAX_APPROX) {
        logger().warn(
            "Enabling the improved max approximator while not using area weighting may lead to incorrect results.");
    }

    m_collision_set_type = type;
}

void NormalCollisions::set_enable_shape_derivatives(
    const bool enable_shape_derivatives)
{
    if (!empty() && enable_shape_derivatives != m_enable_shape_derivatives) {
        logger().warn(
            "Setting enable_shape_derivatives after building collisions. "
            "Re-build collisions for this to have an effect.");
    }

    m_enable_shape_derivatives = enable_shape_derivatives;
}

// ============================================================================

// NOTE: Actually distance squared
double NormalCollisions::compute_minimum_distance(
    const CollisionMesh& mesh, Eigen::ConstRef<Eigen::MatrixXd> vertices) const
{
    assert(vertices.rows() == mesh.num_vertices());

    if (empty()) {
        return std::numeric_limits<double>::infinity();
    }

    const Eigen::MatrixXi& edges = mesh.edges();
    const Eigen::MatrixXi& faces = mesh.faces();

    return tbb::parallel_reduce(
        tbb::blocked_range<size_t>(0, size()),
        std::numeric_limits<double>::infinity(),
        [&](tbb::blocked_range<size_t> r, double partial_min_dist) -> double {
            for (size_t i = r.begin(); i < r.end(); i++) {
                const double dist = (*this)[i].compute_distance(
                    (*this)[i].dof(vertices, edges, faces));
                partial_min_dist = std::min(partial_min_dist, dist);
            }
            return partial_min_dist;
        },
        [](double a, double b) { return std::min(a, b); });
}

// ============================================================================

size_t NormalCollisions::size() const
{
    return vv_collisions.size() + ev_collisions.size() + ee_collisions.size()
        + fv_collisions.size() + pv_collisions.size();
}

bool NormalCollisions::empty() const
{
    return vv_collisions.empty() && ev_collisions.empty()
        && ee_collisions.empty() && fv_collisions.empty()
        && pv_collisions.empty();
}

void NormalCollisions::clear()
{
    vv_collisions.clear();
    ev_collisions.clear();
    ee_collisions.clear();
    fv_collisions.clear();
    pv_collisions.clear();
}

NormalCollision& NormalCollisions::operator[](size_t i)
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

const NormalCollision& NormalCollisions::operator[](size_t i) const
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

bool NormalCollisions::is_vertex_vertex(size_t i) const
{
    return i < vv_collisions.size();
}

bool NormalCollisions::is_edge_vertex(size_t i) const
{
    return i >= vv_collisions.size()
        && i < vv_collisions.size() + ev_collisions.size();
}

bool NormalCollisions::is_edge_edge(size_t i) const
{
    return i >= vv_collisions.size() + ev_collisions.size()
        && i
        < vv_collisions.size() + ev_collisions.size() + ee_collisions.size();
}

bool NormalCollisions::is_face_vertex(size_t i) const
{
    return i
        >= vv_collisions.size() + ev_collisions.size() + ee_collisions.size()
        && i < vv_collisions.size() + ev_collisions.size()
            + ee_collisions.size() + fv_collisions.size();
}

bool NormalCollisions::is_plane_vertex(size_t i) const
{
    return i >= vv_collisions.size() + ev_collisions.size()
            + ee_collisions.size() + fv_collisions.size()
        && i < vv_collisions.size() + ev_collisions.size()
            + ee_collisions.size() + fv_collisions.size()
            + pv_collisions.size();
}

std::string NormalCollisions::to_string(
    const CollisionMesh& mesh, Eigen::ConstRef<Eigen::MatrixXd> vertices) const
{
    std::stringstream ss;
    for (const auto& vv : vv_collisions) {
        ss << "\n"
           << fmt::format(
                  "vv: {} {}, w: {:g}, d: {:g}",
                  std::min(vv.vertex0_id, vv.vertex1_id),
                  std::max(vv.vertex0_id, vv.vertex1_id), vv.weight,
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
        const index_t min_ei = std::min(ee.edge0_id, ee.edge1_id);
        const index_t max_ei = std::max(ee.edge0_id, ee.edge1_id);
        ss << "\n"
           << fmt::format(
                  "ee: {}=({}, {}) {}=({}, {}), w: {:g}, dtype: {}, d: {:g}",
                  min_ei, mesh.edges()(min_ei, 0), mesh.edges()(min_ei, 1),
                  max_ei, mesh.edges()(max_ei, 0), mesh.edges()(max_ei, 1),
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
