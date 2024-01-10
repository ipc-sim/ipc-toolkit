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

template <int dim>
void SmoothCollisions<dim>::build(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const double dhat,
    const double dmin,
    const BroadPhaseMethod broad_phase_method)
{
    VirtualCollisions::build(mesh, vertices, dhat, dmin, broad_phase_method);
}

template <int dim>
std::vector<CandidateType> SmoothCollisions<dim>::get_candidate_types(const int &_dim) const
{
    assert(dim == _dim);
    if constexpr (dim == 2)
    {
        if (quad_type == SurfaceQuadratureType::SinglePoint)
            return {CandidateType::EdgeVertex};
        else
            return {CandidateType::EdgeEdge};
    }
    else
    {
        if (quad_type == SurfaceQuadratureType::SinglePoint)
            return {CandidateType::FaceVertex, CandidateType::EdgeEdge};
        else
        {
            logger().error("3D face-face candidate is not implemented!");
            return {CandidateType::FaceVertex, CandidateType::EdgeEdge};
        }
    }
}

template <int dim>
void SmoothCollisions<dim>::build(
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

    tbb::enumerable_thread_specific<SmoothCollisionsBuilder<dim>> storage;

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), candidates.ev_candidates.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            storage.local().add_edge_vertex_collisions(
                mesh, vertices, candidates.ev_candidates, is_active, r.begin(),
                r.end(), dhat, use_adaptive_eps);
        });

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), candidates.ee_candidates.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            storage.local().add_edge_edge_collisions(
                mesh, vertices, candidates.ee_candidates, is_active, r.begin(),
                r.end(), quad_type, dhat, use_adaptive_eps);
        });

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), candidates.fv_candidates.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            storage.local().add_face_vertex_collisions(
                mesh, vertices, candidates.fv_candidates, is_active, r.begin(),
                r.end());
        });

    // -------------------------------------------------------------------------

    SmoothCollisionsBuilder<dim>::merge(storage, *this);

    // logger().debug(to_string(mesh, vertices));

    for (size_t ci = 0; ci < size(); ci++) {
        Collision& collision = (*this)[ci];
        collision.dmin = dmin;
    }
}

// ============================================================================

template <int dim>
size_t SmoothCollisions<dim>::size() const
{
    return ev_collisions.size() + ee_collisions.size() + fv_collisions.size();
}

template <int dim>
bool SmoothCollisions<dim>::empty() const
{
    return ev_collisions.empty() && ee_collisions.empty() && fv_collisions.empty();
}

template <int dim>
void SmoothCollisions<dim>::clear()
{
    ev_collisions.clear();
    ee_collisions.clear();
    fv_collisions.clear();
}

template <int dim>
Collision& SmoothCollisions<dim>::operator[](size_t i)
{
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
    throw std::out_of_range("Collision index is out of range!");
}

template <int dim>
const Collision& SmoothCollisions<dim>::operator[](size_t i) const
{
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
    throw std::out_of_range("Collision index is out of range!");
}

template <int dim>
bool SmoothCollisions<dim>::is_vertex_vertex(size_t i) const
{
    return false;
}

template <int dim>
bool SmoothCollisions<dim>::is_edge_vertex(size_t i) const
{
    return i < ev_collisions.size();
}

template <int dim>
bool SmoothCollisions<dim>::is_edge_edge(size_t i) const
{
    return i >= ev_collisions.size()
        && i
        < ev_collisions.size() + ee_collisions.size();
}

template <int dim>
bool SmoothCollisions<dim>::is_face_vertex(size_t i) const
{
    return false;
}

template <int dim>
bool SmoothCollisions<dim>::is_plane_vertex(size_t i) const
{
    return i
        >= ev_collisions.size() + ee_collisions.size()
        && i < ev_collisions.size()
            + ee_collisions.size() + fv_collisions.size();
}

template <int dim>
std::string SmoothCollisions<dim>::to_string(
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

template class SmoothCollisions<2>;
template class SmoothCollisions<3>;
} // namespace ipc
