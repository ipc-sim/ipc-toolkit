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
    const std::function<double(const long &)> &vert_dhat,
    const std::function<double(const long &)> &edge_dhat,
    const size_t start_i,
    const size_t end_i)
{
    if constexpr (dim == 2)
    {
        const std::vector<unordered_set<int>>& vertex_edge_adj = mesh.vertex_edge_adjacencies();
        for (size_t i = start_i; i < end_i; i++) {
            const auto& [ei, vi] = candidates[i];

            for (int ej : vertex_edge_adj[vi])
            {
                if (ej != ei)
                {
                    std::array<double, 2> dhats = {{edge_dhat(ej), edge_dhat(ei)}};
                    add_collision<SmoothEdgeEdgeCollision>(std::make_shared<SmoothEdgeEdgeCollision>(ej, ei, mesh, param, dhats, vertices), edge_edge_2_to_id, collisions);
                }
            }
            
            for (int j : {0, 1})
            {
                std::array<double, 2> dhats = {{ vert_dhat(mesh.edges()(ei, j)), vert_dhat(vi) }};
                add_collision<SmoothVertexVertexCollision>(std::make_shared<SmoothVertexVertexCollision>(mesh.edges()(ei, j), vi, mesh, param, dhats, vertices), vert_vert_2_to_id, collisions);
            }
        }
    }
}

template <int dim>
void SmoothCollisionsBuilder<dim>::add_neighbor_edge_collisions(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const ParameterType &param,
        const std::function<double(const long &)> &edge_dhat,
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
                    {
                        std::array<double, 2> dhats = {{edge_dhat(i), edge_dhat(j)}};
                        add_collision<SmoothEdgeEdgeCollision>(std::make_shared<SmoothEdgeEdgeCollision>(i, j, mesh, param, dhats, vertices), edge_edge_2_to_id, collisions);
                    }
    }
}

// ============================================================================

template <int dim>
void SmoothCollisionsBuilder<dim>::add_edge_edge_collisions(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const std::vector<EdgeEdgeCandidate>& candidates,
    const ParameterType &param,
    const std::function<double(const long &)> &edge_dhat,
    const size_t start_i,
    const size_t end_i)
{
    if constexpr (dim == 3)
    {
        // const auto &edges_to_faces_adj = mesh.edges_to_faces();
        for (size_t i = start_i; i < end_i; i++) {
            const auto& [eai, ebi] = candidates[i];

            std::array<double, 2> dhats = {{edge_dhat(eai), edge_dhat(ebi)}};
            add_collision<SmoothEdgeEdge3Collision>(std::make_shared<SmoothEdgeEdge3Collision>(eai, ebi, mesh, param, dhats, vertices), edge_edge_3_to_id, collisions);
        }
    }
}

template <int dim>
void SmoothCollisionsBuilder<dim>::add_face_vertex_collisions(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const std::vector<FaceVertexCandidate>& candidates,
    const ParameterType &param,
    const std::function<double(const long &)> &vert_dhat,
    const std::function<double(const long &)> &face_dhat,
    const size_t start_i,
    const size_t end_i)
{
    if constexpr (dim == 3)
    {
        const auto &vertices_to_faces_adj = mesh.vertices_to_faces();
        for (size_t i = start_i; i < end_i; i++) {
            const auto& [fi, vi] = candidates[i];
            
            for (int fj : vertices_to_faces_adj[vi])
                if (fj != fi)
                {
                    std::array<double, 2> dhats = {{face_dhat(fi), face_dhat(fj)}};
                    add_collision<SmoothFaceFaceCollision>(std::make_shared<SmoothFaceFaceCollision>(fi, fj, mesh, param, dhats, vertices), face_face_to_id, collisions);
                }
            
            for (int lv = 0; lv < 3; lv++)
            {
                const auto &vj = mesh.faces()(fi, lv);
                std::array<double, 2> dhats = {{vert_dhat(vi), vert_dhat(vj)}};
                add_collision<SmoothVertexVertex3Collision>(std::make_shared<SmoothVertexVertex3Collision>(vi, vj, mesh, param, dhats, vertices), vert_vert_3_to_id, collisions);
            }
        }
    }
}

template <int dim>
void SmoothCollisionsBuilder<dim>::add_neighbor_face_collisions(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const ParameterType &param,
        const std::function<double(const long &)> &face_dhat,
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
                    {
                        std::array<double, 2> dhats = {{face_dhat(i), face_dhat(j)}};
                        add_collision<SmoothFaceFaceCollision>(std::make_shared<SmoothFaceFaceCollision>(i, j, mesh, param, dhats, vertices), face_face_to_id, collisions);
                    }
    }
}

template <int dim> template <typename TCollision>
void SmoothCollisionsBuilder<dim>::add_collision(
    const std::shared_ptr<TCollision>& pair,
    unordered_map<unordered_tuple, std::tuple<TCollision, long> >& cc_to_id_,
    std::vector<std::shared_ptr<typename SmoothCollisions<dim>::value_type>>& collisions_)
{
    const auto &cc = *std::dynamic_pointer_cast<TCollision>(pair);
    if (pair->is_active() && cc_to_id_.find(pair->get_hash()) == cc_to_id_.end())
    {
        // New collision, so add it to the end of collisions
        cc_to_id_.emplace(pair->get_hash(), std::tuple<TCollision, long>(cc, collisions_.size()));
        collisions_.push_back(pair);
    }
}

template <int dim>
void SmoothCollisionsBuilder<dim>::merge(
    const tbb::enumerable_thread_specific<SmoothCollisionsBuilder<dim>>& local_storage,
    SmoothCollisions<dim>& merged_collisions)
{
    unordered_map<unordered_tuple, std::tuple<SmoothVertexVertexCollision, long> > vert_vert_2_to_id;
    unordered_map<unordered_tuple, std::tuple<SmoothEdgeEdgeCollision, long> > edge_edge_2_to_id;

    unordered_map<unordered_tuple, std::tuple<SmoothVertexVertex3Collision, long> > vert_vert_3_to_id;
    unordered_map<unordered_tuple, std::tuple<SmoothFaceFaceCollision, long> > face_face_to_id;
    unordered_map<unordered_tuple, std::tuple<SmoothEdgeEdge3Collision, long> > edge_edge_3_to_id;

    // size up the hash items
    size_t total = 0;
    for (const auto& storage : local_storage)
        total += storage.collisions.size();
    
    merged_collisions.collisions.reserve(total);

    // merge
    for (auto& builder : local_storage)
        for (auto& cc : builder.collisions)
        {
            if (auto ee3 = std::dynamic_pointer_cast<SmoothEdgeEdge3Collision>(cc))
            {
                if (edge_edge_3_to_id.find(cc->get_hash()) == edge_edge_3_to_id.end())
                {
                    merged_collisions.collisions.push_back(cc);
                    edge_edge_3_to_id.emplace(cc->get_hash(), std::make_tuple<SmoothEdgeEdge3Collision, long>(std::move(*ee3), merged_collisions.collisions.size()));
                }
            }
            else if (auto vv2 = std::dynamic_pointer_cast<SmoothVertexVertexCollision>(cc))
            {
                if (vert_vert_2_to_id.find(cc->get_hash()) == vert_vert_2_to_id.end())
                {
                    merged_collisions.collisions.push_back(cc);
                    vert_vert_2_to_id.emplace(cc->get_hash(), std::make_tuple<SmoothVertexVertexCollision, long>(std::move(*vv2), merged_collisions.collisions.size()));
                }
            }
            else if (auto vv3 = std::dynamic_pointer_cast<SmoothVertexVertex3Collision>(cc))
            {
                if (vert_vert_3_to_id.find(cc->get_hash()) == vert_vert_3_to_id.end())
                {
                    merged_collisions.collisions.push_back(cc);
                    vert_vert_3_to_id.emplace(cc->get_hash(), std::make_tuple<SmoothVertexVertex3Collision, long>(std::move(*vv3), merged_collisions.collisions.size()));
                }
            }
            else if (auto ee = std::dynamic_pointer_cast<SmoothEdgeEdgeCollision>(cc))
            {
                if (edge_edge_2_to_id.find(cc->get_hash()) == edge_edge_2_to_id.end())
                {
                    merged_collisions.collisions.push_back(cc);
                    edge_edge_2_to_id.emplace(cc->get_hash(), std::make_tuple<SmoothEdgeEdgeCollision, long>(std::move(*ee), merged_collisions.collisions.size()));
                }
            }
            else if (auto ff = std::dynamic_pointer_cast<SmoothFaceFaceCollision>(cc))
            {
                if (face_face_to_id.find(cc->get_hash()) == face_face_to_id.end())
                {
                    merged_collisions.collisions.push_back(cc);
                    face_face_to_id.emplace(cc->get_hash(), std::make_tuple<SmoothFaceFaceCollision, long>(std::move(*ff), merged_collisions.collisions.size()));
                }
            }
            else
                throw std::runtime_error("Invalid collision type!");
        }
}

template class SmoothCollisionsBuilder<2>;
template class SmoothCollisionsBuilder<3>;

} // namespace ipc