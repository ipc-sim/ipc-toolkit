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
    const ParameterType &param,
    const BroadPhaseMethod broad_phase_method)
{
    assert(vertices.rows() == mesh.num_vertices());

    double inflation_radius = sqrt(param.eps) / 2;

    candidates.build(mesh, vertices, inflation_radius, broad_phase_method);
    this->build(candidates, mesh, vertices, param);
}

template <int dim>
void SmoothCollisions<dim>::build(
    const Candidates& candidates_,
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const ParameterType &param)
{
    assert(vertices.rows() == mesh.num_vertices());

    clear();

    tbb::enumerable_thread_specific<SmoothCollisionsBuilder<dim>> storage;
    if constexpr (dim == 2)
    {
        tbb::parallel_for(
            tbb::blocked_range<size_t>(size_t(0), candidates_.ev_candidates.size()),
            [&](const tbb::blocked_range<size_t>& r) {
                storage.local().add_edge_vertex_collisions(
                    mesh, vertices, candidates_.ev_candidates, param, r.begin(),
                    r.end());
            });

        if (use_high_order_quadrature)
            tbb::parallel_for(
                tbb::blocked_range<size_t>(size_t(0), mesh.num_vertices()),
                [&](const tbb::blocked_range<size_t>& r) {
                    storage.local().add_neighbor_edge_collisions(
                        mesh, r.begin(), r.end());
                });
    }
    else
    {
        tbb::parallel_for(
            tbb::blocked_range<size_t>(size_t(0), candidates_.ee_candidates.size()),
            [&](const tbb::blocked_range<size_t>& r) {
                storage.local().add_edge_edge_collisions(
                    mesh, vertices, candidates_.ee_candidates, param, r.begin(),
                    r.end());
            });

        tbb::parallel_for(
            tbb::blocked_range<size_t>(size_t(0), candidates_.fv_candidates.size()),
            [&](const tbb::blocked_range<size_t>& r) {
                storage.local().add_face_vertex_collisions(
                    mesh, vertices, candidates_.fv_candidates, param, r.begin(),
                    r.end());
            });

        if (use_high_order_quadrature)
            tbb::parallel_for(
                tbb::blocked_range<size_t>(size_t(0), mesh.num_vertices()),
                [&](const tbb::blocked_range<size_t>& r) {
                    storage.local().add_neighbor_face_collisions(
                        mesh, vertices, param, r.begin(), r.end());
                });
    }
    SmoothCollisionsBuilder<dim>::merge(storage, *this);
    candidates = candidates_;

    // logger().debug(to_string(mesh, vertices));

    if (use_adaptive_dhat)
        for (size_t ci = 0; ci < size(); ci++) {
            typename SmoothCollisions<dim>::value_type& collision = (*this)[ci];
            collision.set_adaptive_dhat(mesh, sqrt(param.eps));
        }
}

// ============================================================================

template <int dim>
size_t SmoothCollisions<dim>::size() const
{
    return collisions.size();
}

template <int dim>
bool SmoothCollisions<dim>::empty() const
{
    return collisions.empty();
}

template <int dim>
void SmoothCollisions<dim>::clear()
{
    collisions.clear();
}

template <int dim>
typename SmoothCollisions<dim>::value_type& SmoothCollisions<dim>::operator[](size_t i)
{
    if (i < collisions.size()) {
        return *collisions[i];
    }
    throw std::out_of_range("Collision index is out of range!");
}

template <int dim>
const typename SmoothCollisions<dim>::value_type& SmoothCollisions<dim>::operator[](size_t i) const
{
    if (i < collisions.size()) {
        return *collisions[i];
    }
    throw std::out_of_range("Collision index is out of range!");
}

template <int dim>
std::string SmoothCollisions<dim>::to_string(
    const CollisionMesh& mesh, const Eigen::MatrixXd& vertices) const
{
    std::stringstream ss;
    for (const auto& cc : collisions) {
        ss << "\n";
        if (!std::dynamic_pointer_cast<SmoothFaceFaceCollision>(cc))
        {
            ss << fmt::format(
                  "ee: {}=({}, {}) {}=({}, {})",
                  (*cc)[0], mesh.edges()((*cc)[0], 0),
                  mesh.edges()((*cc)[0], 1), (*cc)[1],
                  mesh.edges()((*cc)[1], 0), mesh.edges()((*cc)[1], 1));
        }
        else
        {
            ss << fmt::format(
                  "ff: {}=({}, {}, {}) {}=({}, {}, {})", (*cc)[0],
                  mesh.faces()((*cc)[0], 0), mesh.faces()((*cc)[0], 1),
                  mesh.faces()((*cc)[0], 2), (*cc)[1],
                  mesh.faces()((*cc)[1], 0), mesh.faces()((*cc)[1], 1),
                  mesh.faces()((*cc)[1], 2));
        }
    }
    return ss.str();
}

// NOTE: Actually distance squared
template <int dim>
double SmoothCollisions<dim>::compute_minimum_distance(
    const CollisionMesh& mesh, const Eigen::MatrixXd& vertices) const
{
    return 0.;
    // assert(vertices.rows() == mesh.num_vertices());

    // if (candidates.empty()) {
    //     return std::numeric_limits<double>::infinity();
    // }

    // const Eigen::MatrixXi& edges = mesh.edges();
    // const Eigen::MatrixXi& faces = mesh.faces();

    // tbb::enumerable_thread_specific<double> storage(
    //     std::numeric_limits<double>::infinity());

    // tbb::parallel_for(
    //     tbb::blocked_range<size_t>(0, size()),
    //     [&](tbb::blocked_range<size_t> r) {
    //         double& local_min_dist = storage.local();

    //         for (size_t i = r.begin(); i < r.end(); i++) {
    //             const double dist = (*this)[i].compute_distance(
    //                 (*this)[i].dof(vertices, edges, faces));

    //             if (dist < local_min_dist) {
    //                 local_min_dist = dist;
    //             }
    //         }
    //     });

    // const double min_dist = storage.combine([](double a, double b) { return std::min(a, b); });
    // if (min_dist < 1e-10)
    // {
    //     for (int i = 0; i < size(); ++i)
    //     {
    //         const double dist = (*this)[i].compute_distance(
    //             (*this)[i].dof(vertices, edges, faces));
    //         if (dist <= min_dist * (1 + 1e-12))
    //         {
    //             const std::array<long, 4> idx = (*this)[i].vertex_ids(edges, faces);
    //             std::cout << idx[0] << " " << idx[1] << " " << idx[2] << " " << idx[3] << "\n";
    //             std::cout << vertices.row(idx[0]) << " " << vertices.row(idx[1]) << " " << vertices.row(idx[2]) << "\n";
    //             Eigen::MatrixXd V = vertices;
    //             V.conservativeResize(V.rows(), 3);
    //             V.col(2).setZero();
    //             igl::writePLY("zero-dist.ply", V, faces, edges);
    //             exit(0);
    //         }
    //     }
    // }

    // return min_dist;
}

template class SmoothCollisions<2>;
template class SmoothCollisions<3>;
} // namespace ipc
