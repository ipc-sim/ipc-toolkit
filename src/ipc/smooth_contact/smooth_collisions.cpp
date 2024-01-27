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
void SmoothCollisions<dim>::compute_adaptive_dhat(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices, // set to zero for rest pose
    const ParameterType &param,
    const BroadPhaseMethod broad_phase_method)
{
    assert(vertices.rows() == mesh.num_vertices());

    const double dhat = param.dhat;
    double inflation_radius = dhat / 2;

    // Candidates candidates;
    candidates.build(mesh, vertices, inflation_radius, broad_phase_method);
    this->build(candidates, mesh, vertices, param, false /*disable adaptive dhat to compute true pairs*/);

    vert_adaptive_dhat.resize(mesh.num_vertices());
    vert_adaptive_dhat.setConstant(dhat);
    edge_adaptive_dhat.resize(mesh.num_edges());
    edge_adaptive_dhat.setConstant(dhat);
    face_adaptive_dhat.resize(mesh.num_faces());
    face_adaptive_dhat.setConstant(dhat);

    for (auto cc : collisions)
    {
        const double dist = param.get_adaptive_dhat_ratio() * sqrt(cc->compute_distance(cc->dof(vertices, mesh.edges(), mesh.faces())));
        if (std::dynamic_pointer_cast<SmoothEdgeEdge3Collision>(cc))
        {
            edge_adaptive_dhat((*cc)[0]) = std::min(edge_adaptive_dhat((*cc)[0]), dist);
            edge_adaptive_dhat((*cc)[1]) = std::min(edge_adaptive_dhat((*cc)[1]), dist);
        }
        else if (std::dynamic_pointer_cast<SmoothEdgeVertexCollision>(cc) || std::dynamic_pointer_cast<SmoothEdgeVertex3Collision>(cc))
        {
            edge_adaptive_dhat((*cc)[0]) = std::min(edge_adaptive_dhat((*cc)[0]), dist);
            vert_adaptive_dhat((*cc)[1]) = std::min(vert_adaptive_dhat((*cc)[1]), dist);
        }
        else if (std::dynamic_pointer_cast<SmoothFaceVertexCollision>(cc))
        {
            face_adaptive_dhat((*cc)[0]) = std::min(face_adaptive_dhat((*cc)[0]), dist);
            vert_adaptive_dhat((*cc)[1]) = std::min(vert_adaptive_dhat((*cc)[1]), dist);
        }
        else if (std::dynamic_pointer_cast<SmoothVertexVertexCollision>(cc) || std::dynamic_pointer_cast<SmoothVertexVertex3Collision>(cc))
        {
            vert_adaptive_dhat((*cc)[0]) = std::min(vert_adaptive_dhat((*cc)[0]), dist);
            vert_adaptive_dhat((*cc)[1]) = std::min(vert_adaptive_dhat((*cc)[1]), dist);
        }
        else
            throw std::runtime_error("Invalid collision type!");
    }

    for (int f = 0; f < mesh.num_faces(); f++)
        for (int lv = 0; lv < 3; lv++)
        {
            face_adaptive_dhat(f) = std::min(face_adaptive_dhat(f), vert_adaptive_dhat(mesh.faces()(f, lv)));
            face_adaptive_dhat(f) = std::min(face_adaptive_dhat(f), edge_adaptive_dhat(mesh.faces_to_edges()(f, lv)));
        }

    for (int e = 0; e < mesh.num_edges(); e++)
        for (int lv = 0; lv < 2; lv++)
            edge_adaptive_dhat(e) = std::min(edge_adaptive_dhat(e), vert_adaptive_dhat(mesh.edges()(e, lv)));

    logger().debug("vert dhat min {:.2e}, max {:.2e}", vert_adaptive_dhat.minCoeff(), vert_adaptive_dhat.maxCoeff());
    logger().debug("edge dhat min {:.2e}, max {:.2e}", edge_adaptive_dhat.minCoeff(), edge_adaptive_dhat.maxCoeff());
    if (mesh.dim() == 3)
        logger().debug("face dhat min {:.2e}, max {:.2e}", face_adaptive_dhat.minCoeff(), face_adaptive_dhat.maxCoeff());
}

template <int dim>
void SmoothCollisions<dim>::build(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const ParameterType &param,
    const bool use_adaptive_dhat,
    const BroadPhaseMethod broad_phase_method)
{
    assert(vertices.rows() == mesh.num_vertices());

    double inflation_radius = param.dhat / 2;

    // Candidates candidates;
    candidates.build(mesh, vertices, inflation_radius, broad_phase_method);
    this->build(candidates, mesh, vertices, param, use_adaptive_dhat);
}

template <int dim>
void SmoothCollisions<dim>::build(
    const Candidates& candidates_,
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const ParameterType &param,
    const bool use_adaptive_dhat)
{
    assert(vertices.rows() == mesh.num_vertices());

    clear();

    const double dhat = param.dhat;
    if (!use_adaptive_dhat)
    {
        vert_adaptive_dhat.resize(1);
        vert_adaptive_dhat(0) = dhat;
        edge_adaptive_dhat.resize(1);
        edge_adaptive_dhat(0) = dhat;
        face_adaptive_dhat.resize(1);
        face_adaptive_dhat(0) = dhat;
    }

    auto vert_dhat = [&](const long &v_id) {
        return this->get_vert_dhat(v_id);
    };
    auto edge_dhat = [&](const long &e_id) {
        return this->get_edge_dhat(e_id);
    };
    auto face_dhat = [&](const long &f_id) {
        return this->get_face_dhat(f_id);
    };

    tbb::enumerable_thread_specific<SmoothCollisionsBuilder<dim>> storage;
    if constexpr (dim == 2)
    {
        tbb::parallel_for(
            tbb::blocked_range<size_t>(size_t(0), candidates_.ev_candidates.size()),
            [&](const tbb::blocked_range<size_t>& r) {
                storage.local().add_edge_vertex_collisions(
                    mesh, vertices, candidates_.ev_candidates, param, vert_dhat, edge_dhat, r.begin(),
                    r.end());
            });

        // if (use_high_order_quadrature)
        //     tbb::parallel_for(
        //         tbb::blocked_range<size_t>(size_t(0), mesh.num_vertices()),
        //         [&](const tbb::blocked_range<size_t>& r) {
        //             storage.local().add_neighbor_edge_collisions(
        //                 mesh, vertices, param, edge_dhat, r.begin(), r.end());
        //         });
    }
    else
    {
        tbb::parallel_for(
            tbb::blocked_range<size_t>(size_t(0), candidates_.ee_candidates.size()),
            [&](const tbb::blocked_range<size_t>& r) {
                storage.local().add_edge_edge_collisions(
                    mesh, vertices, candidates_.ee_candidates, param, vert_dhat, edge_dhat, r.begin(),
                    r.end());
            });

        tbb::parallel_for(
            tbb::blocked_range<size_t>(size_t(0), candidates_.fv_candidates.size()),
            [&](const tbb::blocked_range<size_t>& r) {
                storage.local().add_face_vertex_collisions(
                    mesh, vertices, candidates_.fv_candidates, param, vert_dhat, edge_dhat, face_dhat, r.begin(),
                    r.end());
            });

        // if (use_high_order_quadrature)
        //     tbb::parallel_for(
        //         tbb::blocked_range<size_t>(size_t(0), mesh.num_vertices()),
        //         [&](const tbb::blocked_range<size_t>& r) {
        //             storage.local().add_neighbor_face_collisions(
        //                 mesh, vertices, param, face_dhat, r.begin(), r.end());
        //         });
    }
    SmoothCollisionsBuilder<dim>::merge(storage, *this);
    candidates = candidates_;

    // logger().debug(to_string(mesh, vertices, param));
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
    const CollisionMesh& mesh, const Eigen::MatrixXd& vertices, const ParameterType &params) const
{
    std::stringstream ss;
    for (const auto& cc : collisions) {
        ss << "\n";
        {
            ss << fmt::format(
                  "[{}]: ({} {}) dist {} potential {} grad {}", cc->name(), (*cc)[0], (*cc)[1],
                    cc->compute_distance(cc->dof(vertices, mesh.edges(), mesh.faces())),
                    (*cc)(cc->dof(vertices, mesh.edges(), mesh.faces()), params),
                    (*cc).gradient(cc->dof(vertices, mesh.edges(), mesh.faces()), params).norm()
                  );
        }
    }
    return ss.str();
}

// NOTE: Actually distance squared
template <int dim>
double SmoothCollisions<dim>::compute_minimum_distance(
    const CollisionMesh& mesh, const Eigen::MatrixXd& vertices) const
{
    assert(vertices.rows() == mesh.num_vertices());

    if (candidates.empty()) {
        return std::numeric_limits<double>::infinity();
    }

    const Eigen::MatrixXi& edges = mesh.edges();
    const Eigen::MatrixXi& faces = mesh.faces();

    tbb::enumerable_thread_specific<double> storage(
        std::numeric_limits<double>::infinity());

    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, candidates.size()),
        [&](tbb::blocked_range<size_t> r) {
            double& local_min_dist = storage.local();

            for (size_t i = r.begin(); i < r.end(); i++) {
                const double dist = candidates[i].compute_distance(
                    candidates[i].dof(vertices, edges, faces));

                if (dist < local_min_dist) {
                    local_min_dist = dist;
                }
            }
        });

    return storage.combine([](double a, double b) { return std::min(a, b); });
}

template class SmoothCollisions<2>;
template class SmoothCollisions<3>;
} // namespace ipc
