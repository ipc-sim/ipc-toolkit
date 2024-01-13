#include "smooth_collisions_builder.hpp"

#include <ipc/distance/distance_type.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/edge_edge_mollifier.hpp>
#include <ipc/distance/point_triangle.hpp>

namespace ipc {

template <int dim>
void SmoothCollisionsBuilder<dim>::add_edge_vertex_collisions(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const std::vector<EdgeVertexCandidate>& candidates,
    const std::function<bool(double)>& is_active,
    const size_t start_i,
    const size_t end_i,
    const double dhat)
{
    if constexpr (dim == 2)
    {
        const std::vector<unordered_set<int>>& vertex_edge_adj = mesh.vertex_edge_adjacencies();
        for (size_t i = start_i; i < end_i; i++) {
            const auto& [ei, vi] = candidates[i];
            const auto [v, e0, e1, _] =
                candidates[i].vertices(vertices, mesh.edges(), mesh.faces());
            const PointEdgeDistanceType dtype = point_edge_distance_type(v, e0, e1);
            const double distance_sqr = point_edge_distance(v, e0, e1, dtype);

            if (!is_active(distance_sqr))
                continue;
            
            for (int ej : vertex_edge_adj[vi])
                if (ej != ei)
                    add_collision<SmoothEdgeEdgeCollision>(mesh, unordered_tuple(ej, ei), cc_to_id, collisions);
        }
    }
}

template <int dim>
void SmoothCollisionsBuilder<dim>::add_neighbor_edge_collisions(
        const CollisionMesh& mesh,
        const size_t start_i,
        const size_t end_i)
{
    if constexpr (dim == 2)
    {
        const std::vector<unordered_set<int>>& vertex_edge_adj = mesh.vertex_edge_adjacencies();
        for (size_t v = start_i; v < end_i; v++)
            for (int i : vertex_edge_adj[v])
                for (int j : vertex_edge_adj[v])
                    if (j > i)
                        add_collision<SmoothEdgeEdgeCollision>(mesh, unordered_tuple(i, j), cc_to_id, collisions);
    }
}

// ============================================================================

template <int dim>
void SmoothCollisionsBuilder<dim>::add_edge_edge_collisions(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const std::vector<EdgeEdgeCandidate>& candidates,
    const std::function<bool(double)>& is_active,
    const size_t start_i,
    const size_t end_i,
    const double dhat)
{
    if constexpr (dim == 3)
    {
        // const auto &edges_to_faces_adj = mesh.edges_to_faces();
        for (size_t i = start_i; i < end_i; i++) {
            const auto& [eai, ebi] = candidates[i];

            const auto [ea0, ea1, eb0, eb1] =
                candidates[i].vertices(vertices, mesh.edges(), mesh.faces());

            const double distance_sqr = edge_edge_distance(ea0, ea1, eb0, eb1, EdgeEdgeDistanceType::AUTO);

            if (!is_active(distance_sqr))
                continue;

            add_collision<SmoothEdgeEdge3Collision>(mesh, unordered_tuple(eai, ebi), ee_to_id, collisions);
        }
    }
}

template <int dim>
void SmoothCollisionsBuilder<dim>::add_face_vertex_collisions(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const std::vector<FaceVertexCandidate>& candidates,
    const std::function<bool(double)>& is_active,
    const size_t start_i,
    const size_t end_i)
{
    if constexpr (dim == 3)
    {
        const auto &vertices_to_faces_adj = mesh.vertices_to_faces();
        for (size_t i = start_i; i < end_i; i++) {
            const auto& [fi, vi] = candidates[i];

            const auto [v, f0, f1, f2] =
                candidates[i].vertices(vertices, mesh.edges(), mesh.faces());

            // Compute distance type
            const PointTriangleDistanceType dtype =
                point_triangle_distance_type(v, f0, f1, f2);
            const double distance_sqr =
                point_triangle_distance(v, f0, f1, f2, dtype);

            if (!is_active(distance_sqr))
                continue;
            
            for (int fj : vertices_to_faces_adj[vi])
                if (fj != fi)
                    add_collision<SmoothFaceFaceCollision>(mesh, unordered_tuple(fi, fj), cc_to_id, collisions);
        }
    }
}

template <int dim>
void SmoothCollisionsBuilder<dim>::add_neighbor_face_collisions(
        const CollisionMesh& mesh,
        const size_t start_i,
        const size_t end_i)
{
    if constexpr (dim == 3)
    {
        const auto &vertices_to_faces_adj = mesh.vertices_to_faces();
        for (size_t v = start_i; v < end_i; v++)
            for (int i : vertices_to_faces_adj[v])
                for (int j : vertices_to_faces_adj[v])
                    if (j > i)
                        add_collision<SmoothFaceFaceCollision>(mesh, unordered_tuple(i, j), cc_to_id, collisions);
    }
}

template <int dim> template <typename TCollision>
void SmoothCollisionsBuilder<dim>::add_collision(
    const CollisionMesh &mesh,
    const unordered_tuple& pair,
    unordered_map<unordered_tuple, long>& cc_to_id_,
    std::vector<std::shared_ptr<typename SmoothCollisions<dim>::value_type>>& collisions_)
{
    if (cc_to_id_.find(pair) == cc_to_id_.end())
    {
        // New collision, so add it to the end of collisions
        cc_to_id_.emplace(pair, collisions_.size());
        collisions_.push_back(std::make_shared<TCollision>(pair[0], pair[1], mesh));
    }
}

template <int dim>
void SmoothCollisionsBuilder<dim>::merge(
    const tbb::enumerable_thread_specific<SmoothCollisionsBuilder<dim>>& local_storage,
    SmoothCollisions<dim>& merged_collisions)
{
    unordered_map<unordered_tuple, long> cc_to_id, ee_to_id;

    // size up the hash items
    size_t total = 0;
    for (const auto& storage : local_storage)
        total += storage.collisions.size();
    
    merged_collisions.collisions.reserve(total);

    // merge
    for (const auto& builder : local_storage)
        for (const auto& cc : builder.collisions)
        {
            unordered_tuple pair(cc->primitive0, cc->primitive1);
            if (std::dynamic_pointer_cast<SmoothEdgeEdge3Collision>(cc))
            {
                if (ee_to_id.find(pair) == ee_to_id.end())
                {
                    // New collision, so add it to the end of collisions
                    ee_to_id.emplace(pair, merged_collisions.collisions.size());
                    merged_collisions.collisions.push_back(cc);
                }
            }
            else if (cc_to_id.find(pair) == cc_to_id.end())
            {
                // New collision, so add it to the end of collisions
                cc_to_id.emplace(pair, merged_collisions.collisions.size());
                merged_collisions.collisions.push_back(cc);
            }
        }
}

template class SmoothCollisionsBuilder<2>;
template class SmoothCollisionsBuilder<3>;

} // namespace ipc