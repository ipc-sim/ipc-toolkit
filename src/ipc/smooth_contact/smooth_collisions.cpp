#include "smooth_collisions.hpp"

#include "smooth_collisions_builder.hpp"
#include <ipc/distance/point_line.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/distance/edge_edge.hpp>
#include <ipc/utils/local_to_global.hpp>

#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>
#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>

#include <stdexcept> // std::out_of_range

namespace ipc {

void SmoothCollisions::build(
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

void SmoothCollisions::build(
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

    tbb::enumerable_thread_specific<SmoothCollisionsBuilder> storage;

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

    // -------------------------------------------------------------------------

    SmoothCollisionsBuilder::merge(storage, *this);

    // logger().debug(to_string(mesh, vertices));

    for (size_t ci = 0; ci < size(); ci++) {
        Collision& collision = (*this)[ci];
        collision.dmin = dmin;
    }
}

// ============================================================================

// NOTE: Actually distance squared
double SmoothCollisions::compute_minimum_distance(
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

size_t SmoothCollisions::size() const
{
    return ev_collisions.size();
}

bool SmoothCollisions::empty() const
{
    return ev_collisions.empty();
}

void SmoothCollisions::clear()
{
    ev_collisions.clear();
}

Collision& SmoothCollisions::operator[](size_t i)
{
    if (i < ev_collisions.size()) {
        return ev_collisions[i];
    }
    i -= ev_collisions.size();
    throw std::out_of_range("Collision index is out of range!");
}

const Collision& SmoothCollisions::operator[](size_t i) const
{
    if (i < ev_collisions.size()) {
        return ev_collisions[i];
    }
    i -= ev_collisions.size();
    throw std::out_of_range("Collision index is out of range!");
}

bool SmoothCollisions::is_vertex_vertex(size_t i) const
{
    return false;
}

bool SmoothCollisions::is_edge_vertex(size_t i) const
{
    return i < ev_collisions.size();
}

bool SmoothCollisions::is_edge_edge(size_t i) const
{
    return false;
}

bool SmoothCollisions::is_face_vertex(size_t i) const
{
    return false;
}

bool SmoothCollisions::is_plane_vertex(size_t i) const
{
    return false;
}

std::string SmoothCollisions::to_string(
    const CollisionMesh& mesh, const Eigen::MatrixXd& vertices) const
{
    std::stringstream ss;
    for (const auto& ev : ev_collisions) {
        ss << "\n"
           << fmt::format(
                  "ev: {}=({}, {}) {}, w: {:g}, d: {:g}", ev.edge_id,
                  mesh.edges()(ev.edge_id, 0), mesh.edges()(ev.edge_id, 1),
                  ev.vertex_id, ev.weight,
                  ev.compute_distance(
                      ev.dof(vertices, mesh.edges(), mesh.faces())));
    }
    return ss.str();
}

} // namespace ipc
