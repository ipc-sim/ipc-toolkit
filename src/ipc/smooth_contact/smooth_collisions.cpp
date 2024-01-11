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

template <int dim, class TCollision>
void SmoothCollisions<dim, TCollision>::build(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const double dhat,
    const double dmin,
    const BroadPhaseMethod broad_phase_method)
{
    VirtualCollisions<max_vert>::build(mesh, vertices, dhat, dmin, broad_phase_method);
}

template <int dim, class TCollision>
void SmoothCollisions<dim, TCollision>::build(
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

    tbb::enumerable_thread_specific<SmoothCollisionsBuilder<dim, TCollision>> storage;
    if constexpr (dim == 2)
    {
        tbb::parallel_for(
            tbb::blocked_range<size_t>(size_t(0), candidates.ev_candidates.size()),
            [&](const tbb::blocked_range<size_t>& r) {
                storage.local().add_edge_vertex_collisions(
                    mesh, vertices, candidates.ev_candidates, is_active, r.begin(),
                    r.end(), dhat, use_adaptive_eps);
            });

        tbb::parallel_for(
            tbb::blocked_range<size_t>(size_t(0), mesh.num_edges()),
            [&](const tbb::blocked_range<size_t>& r) {
                storage.local().add_neighbor_edge_collisions(
                    mesh, r.begin(), r.end(), use_adaptive_eps);
            });
    }
    else
    {
        tbb::parallel_for(
            tbb::blocked_range<size_t>(size_t(0), candidates.ee_candidates.size()),
            [&](const tbb::blocked_range<size_t>& r) {
                storage.local().add_edge_edge_collisions(
                    mesh, vertices, candidates.ee_candidates, is_active, r.begin(),
                    r.end(), dhat, use_adaptive_eps);
            });

        tbb::parallel_for(
            tbb::blocked_range<size_t>(size_t(0), candidates.fv_candidates.size()),
            [&](const tbb::blocked_range<size_t>& r) {
                storage.local().add_face_vertex_collisions(
                    mesh, vertices, candidates.fv_candidates, is_active, r.begin(),
                    r.end());
            });
    }
    SmoothCollisionsBuilder<dim, TCollision>::merge(storage, *this);

    // logger().debug(to_string(mesh, vertices));

    for (size_t ci = 0; ci < size(); ci++) {
        typename SmoothCollisions<dim, TCollision>::value_type& collision = (*this)[ci];
        collision.dmin = dmin;
    }
}

// ============================================================================

template <int dim, class TCollision>
size_t SmoothCollisions<dim, TCollision>::size() const
{
    return collisions.size();
}

template <int dim, class TCollision>
bool SmoothCollisions<dim, TCollision>::empty() const
{
    return collisions.empty();
}

template <int dim, class TCollision>
void SmoothCollisions<dim, TCollision>::clear()
{
    collisions.clear();
}

template <int dim, class TCollision>
typename SmoothCollisions<dim, TCollision>::value_type& SmoothCollisions<dim, TCollision>::operator[](size_t i)
{
    if (i < collisions.size()) {
        return *collisions[i];
    }
    throw std::out_of_range("Collision index is out of range!");
}

template <int dim, class TCollision>
const typename SmoothCollisions<dim, TCollision>::value_type& SmoothCollisions<dim, TCollision>::operator[](size_t i) const
{
    if (i < collisions.size()) {
        return *collisions[i];
    }
    throw std::out_of_range("Collision index is out of range!");
}

template <int dim, class TCollision>
std::string SmoothCollisions<dim, TCollision>::to_string(
    const CollisionMesh& mesh, const Eigen::MatrixXd& vertices) const
{
    std::stringstream ss;
    for (const auto& cc : collisions) {
        ss << "\n";
        if constexpr (dim == 2)
        {
            ss << fmt::format(
                  "ee: {}=({}, {}) {}=({}, {}), w: {:g}, dtype: {}, d: {:g}",
                  cc->edge0_id, mesh.edges()(cc->edge0_id, 0),
                  mesh.edges()(cc->edge0_id, 1), cc->edge1_id,
                  mesh.edges()(cc->edge1_id, 0), mesh.edges()(cc->edge1_id, 1),
                  cc->weight, int(cc->dtype),
                  cc->compute_distance(
                      cc->dof(vertices, mesh.edges(), mesh.faces())));
        }
        else
        {
            ss << fmt::format(
                  "ff: {}=({}, {}, {}) {}=({}, {}, {}), w: {:g}, d: {:g}", cc->face0_id,
                  mesh.faces()(cc->face0_id, 0), mesh.faces()(cc->face0_id, 1),
                  mesh.faces()(cc->face0_id, 2), cc->face1_id,
                  mesh.faces()(cc->face1_id, 0), mesh.faces()(cc->face1_id, 1),
                  mesh.faces()(cc->face1_id, 2), cc->weight,
                  cc->compute_distance(
                      cc->dof(vertices, mesh.edges(), mesh.faces())));
        }
    }
    return ss.str();
}

template class SmoothCollisions<2, SmoothEdgeEdgeCollision<2>>;
template class SmoothCollisions<3, SmoothFaceFaceCollision>;
} // namespace ipc
