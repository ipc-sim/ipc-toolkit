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
    const ParameterType &param,
    const size_t start_i,
    const size_t end_i)
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

            if (distance_sqr >= param.eps)
                continue;
            
            for (int ej : vertex_edge_adj[vi])
                if (ej != ei)
                    add_collision<SmoothEdgeEdgeCollision>(std::make_shared<SmoothEdgeEdgeCollision>(ej, ei, mesh, param, vertices), edge_edge_2_to_id, collisions);
        }
    }
}

template <int dim>
void SmoothCollisionsBuilder<dim>::add_neighbor_edge_collisions(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const ParameterType &param,
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
                        add_collision<SmoothEdgeEdgeCollision>(std::make_shared<SmoothEdgeEdgeCollision>(i, j, mesh, param, vertices), edge_edge_2_to_id, collisions);
    }
}

// ============================================================================

template <int dim>
void SmoothCollisionsBuilder<dim>::add_edge_edge_collisions(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const std::vector<EdgeEdgeCandidate>& candidates,
    const ParameterType &param,
    const size_t start_i,
    const size_t end_i)
{
    if constexpr (dim == 3)
    {
        // const auto &edges_to_faces_adj = mesh.edges_to_faces();
        for (size_t i = start_i; i < end_i; i++) {
            const auto& [eai, ebi] = candidates[i];

            const auto [ea0, ea1, eb0, eb1] =
                candidates[i].vertices(vertices, mesh.edges(), mesh.faces());

            const double distance_sqr = edge_edge_distance(ea0, ea1, eb0, eb1, EdgeEdgeDistanceType::AUTO);

            if (distance_sqr >= param.eps)
                continue;

            add_collision<SmoothEdgeEdge3Collision>(std::make_shared<SmoothEdgeEdge3Collision>(eai, ebi, mesh, param, vertices), edge_edge_3_to_id, collisions);
        }
    }
}

template <int dim>
void SmoothCollisionsBuilder<dim>::add_face_vertex_collisions(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const std::vector<FaceVertexCandidate>& candidates,
    const ParameterType &param,
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

            if (distance_sqr >= param.eps)
                continue;
            
            for (int fj : vertices_to_faces_adj[vi])
                if (fj != fi)
                    add_collision<SmoothFaceFaceCollision>(std::make_shared<SmoothFaceFaceCollision>(fi, fj, mesh, param, vertices), face_face_to_id, collisions);
        }
    }
}

template <int dim>
void SmoothCollisionsBuilder<dim>::add_neighbor_face_collisions(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const ParameterType &param,
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
                        add_collision<SmoothFaceFaceCollision>(std::make_shared<SmoothFaceFaceCollision>(i, j, mesh, param, vertices), face_face_to_id, collisions);
    }
}

template <int dim> template <typename TCollision>
void SmoothCollisionsBuilder<dim>::add_collision(
    const std::shared_ptr<TCollision>& pair,
    unordered_map<TCollision, long>& cc_to_id_,
    std::vector<std::shared_ptr<typename SmoothCollisions<dim>::value_type>>& collisions_)
{
    const auto &cc = *std::dynamic_pointer_cast<TCollision>(pair);
    if (pair->is_active() && cc_to_id_.find(cc) == cc_to_id_.end())
    {
        // New collision, so add it to the end of collisions
        cc_to_id_.emplace(*pair, collisions_.size());
        collisions_.push_back(pair);
    }
}

template <int dim>
void SmoothCollisionsBuilder<dim>::merge(
    const tbb::enumerable_thread_specific<SmoothCollisionsBuilder<dim>>& local_storage,
    SmoothCollisions<dim>& merged_collisions)
{
    unordered_map<SmoothEdgeEdgeCollision, long> edge_edge_2_to_id;
    unordered_map<SmoothFaceFaceCollision, long> face_face_to_id;
    unordered_map<SmoothEdgeEdge3Collision, long> edge_edge_3_to_id;

    // size up the hash items
    size_t total = 0;
    for (const auto& storage : local_storage)
        total += storage.collisions.size();
    
    merged_collisions.collisions.reserve(total);

    // merge
    for (const auto& builder : local_storage)
        for (const auto& cc : builder.collisions)
        {
            if (auto ee3 = std::dynamic_pointer_cast<SmoothEdgeEdge3Collision>(cc))
            {
                if (edge_edge_3_to_id.find(*ee3) == edge_edge_3_to_id.end())
                {
                    // New collision, so add it to the end of collisions
                    edge_edge_3_to_id.emplace(*ee3, merged_collisions.collisions.size());
                    merged_collisions.collisions.push_back(cc);
                }
            }
            else if (auto ee = std::dynamic_pointer_cast<SmoothEdgeEdgeCollision>(cc))
            {
                if (edge_edge_2_to_id.find(*ee) == edge_edge_2_to_id.end())
                {
                    // New collision, so add it to the end of collisions
                    edge_edge_2_to_id.emplace(*ee, merged_collisions.collisions.size());
                    merged_collisions.collisions.push_back(cc);
                }
            }
            else if (auto ff = std::dynamic_pointer_cast<SmoothFaceFaceCollision>(cc))
            {
                if (face_face_to_id.find(*ff) == face_face_to_id.end())
                {
                    // New collision, so add it to the end of collisions
                    face_face_to_id.emplace(*ff, merged_collisions.collisions.size());
                    merged_collisions.collisions.push_back(cc);
                }
            }
            else
                throw std::runtime_error("Invalid collision type!");
        }
}

template class SmoothCollisionsBuilder<2>;
template class SmoothCollisionsBuilder<3>;

} // namespace ipc