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
        // const std::vector<unordered_set<int>>& vertex_edge_adj = mesh.vertex_edge_adjacencies();
        for (size_t i = start_i; i < end_i; i++) {
            const auto& [ei, vi] = candidates[i];

            // for (int ej : vertex_edge_adj[vi])
            {
                // if (ej != ei)
                {
                    std::array<double, 2> dhats = {{edge_dhat(ei), vert_dhat(vi)}};
                    add_collision<SmoothEdgeVertexCollision>(std::make_shared<SmoothEdgeVertexCollision>(ei, vi, mesh, param, dhats, vertices), vert_edge_2_to_id, collisions);
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

// template <int dim>
// void SmoothCollisionsBuilder<dim>::add_neighbor_edge_collisions(
//         const CollisionMesh& mesh,
//         const Eigen::MatrixXd& vertices,
//         const ParameterType &param,
//         const std::function<double(const long &)> &edge_dhat,
//         const size_t start_i,
//         const size_t end_i)
// {
//     if constexpr (dim == 2)
//     {
//         const std::vector<unordered_set<int>>& vertex_edge_adj = mesh.vertex_edge_adjacencies();
//         for (size_t v = start_i; v < end_i; v++)
//             for (int i : vertex_edge_adj[v])
//                 for (int j : vertex_edge_adj[v])
//                     if (j > i)
//                     {
//                         std::array<double, 2> dhats = {{edge_dhat(i), edge_dhat(j)}};
//                         add_collision<SmoothEdgeVertexCollision>(std::make_shared<SmoothEdgeVertexCollision>(i, j, mesh, param, dhats, vertices), vert_edge_2_to_id, collisions);
//                     }
//     }
// }

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
        for (size_t i = start_i; i < end_i; i++) {
            const auto& [eai, ebi] = candidates[i];

            // {
            //     const auto [ea0i, ea1i, eb0i, eb1i] =
            //         candidates[i].vertex_ids(mesh.edges(), mesh.faces());

            //     const auto [ea0, ea1, eb0, eb1] =
            //         candidates[i].vertices(vertices, mesh.edges(), mesh.faces());

            //     const EdgeEdgeDistanceType actual_dtype =
            //         edge_edge_distance_type(ea0, ea1, eb0, eb1);

            //     const double distance_sqr =
            //         edge_edge_distance(ea0, ea1, eb0, eb1, actual_dtype);
                
            //     if (distance_sqr < 1e-10)
            //         logger().error("edge {} {} dist {}", eai, ebi, sqrt(distance_sqr));
            // }

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
    const std::function<double(const long &)> &edge_dhat,
    const std::function<double(const long &)> &face_dhat,
    const size_t start_i,
    const size_t end_i)
{
    if constexpr (dim == 3)
    {
        // const auto &vertices_to_faces_adj = mesh.vertices_to_faces();
        for (size_t i = start_i; i < end_i; i++) {
            const auto& [fi, vi] = candidates[i];
            
            {
                std::array<double, 2> dhats = {{face_dhat(fi), vert_dhat(vi)}};

                Eigen::Vector3d v = vertices.row(vi);
                Eigen::Vector3d f0 = vertices.row(mesh.faces()(fi, 0));
                Eigen::Vector3d f1 = vertices.row(mesh.faces()(fi, 1));
                Eigen::Vector3d f2 = vertices.row(mesh.faces()(fi, 2));
                const PointTriangleDistanceType dtype =
                    point_triangle_distance_type(v, f0, f1, f2);
                const double distance_sqr =
                    point_triangle_distance(v, f0, f1, f2, dtype);

                if (distance_sqr >= pow(std::min(dhats[0], dhats[1]), 2))
                    continue;

                add_collision<SmoothFaceVertexCollision>(std::make_shared<SmoothFaceVertexCollision>(fi, vi, mesh, param, dhats, vertices), face_vert_to_id, collisions);
            }
            
            for (int lv = 0; lv < 3; lv++)
            {
                const auto &vj = mesh.faces()(fi, lv);
                std::array<double, 2> dhats = {{vert_dhat(vi), vert_dhat(vj)}};
                if ((vertices.row(vi) - vertices.row(vj)).norm() >= std::min(dhats[0], dhats[1]))
                    continue;
                add_collision<SmoothVertexVertex3Collision>(std::make_shared<SmoothVertexVertex3Collision>(vi, vj, mesh, param, dhats, vertices), vert_vert_3_to_id, collisions);
            }

            for (int le = 0; le < 3; le++)
            {
                const auto &eid = mesh.faces_to_edges()(fi, le);
                std::array<double, 2> dhats = {{edge_dhat(eid), vert_dhat(vi)}};

                const PointEdgeDistanceType dtype = point_edge_distance_type(vertices.row(vi), vertices.row(mesh.edges()(eid, 0)), vertices.row(mesh.edges()(eid, 1)));
                const double distance_sqr = point_edge_distance(vertices.row(vi), vertices.row(mesh.edges()(eid, 0)), vertices.row(mesh.edges()(eid, 1)), dtype);

                if (distance_sqr >= pow(std::min(dhats[0], dhats[1]), 2))
                    continue;

                add_collision<SmoothEdgeVertex3Collision>(std::make_shared<SmoothEdgeVertex3Collision>(eid, vi, mesh, param, dhats, vertices), edge_vert_3_to_id, collisions);
            }
        }
    }
}

// template <int dim>
// void SmoothCollisionsBuilder<dim>::add_neighbor_face_collisions(
//         const CollisionMesh& mesh,
//         const Eigen::MatrixXd& vertices,
//         const ParameterType &param,
//         const std::function<double(const long &)> &face_dhat,
//         const size_t start_i,
//         const size_t end_i)
// {
//     if constexpr (dim == 3)
//     {
//         const auto &vertices_to_faces_adj = mesh.vertices_to_faces();
//         for (size_t v = start_i; v < end_i; v++)
//             for (int i : vertices_to_faces_adj[v])
//                 for (int j : vertices_to_faces_adj[v])
//                     if (j > i)
//                     {
//                         std::array<double, 2> dhats = {{face_dhat(i), face_dhat(j)}};
//                         add_collision<SmoothFaceVertexCollision>(std::make_shared<SmoothFaceVertexCollision>(i, j, mesh, param, dhats, vertices), face_face_to_id, collisions);
//                     }
//     }
// }

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
    unordered_map<unordered_tuple, std::tuple<SmoothEdgeVertexCollision, long> > vert_edge_2_to_id;
    
    unordered_map<unordered_tuple, std::tuple<SmoothFaceVertexCollision, long> > face_vert_to_id;
    unordered_map<unordered_tuple, std::tuple<SmoothVertexVertex3Collision, long> > vert_vert_3_to_id;
    unordered_map<unordered_tuple, std::tuple<SmoothEdgeVertex3Collision, long> > edge_vert_3_to_id;
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
            else if (auto ev = std::dynamic_pointer_cast<SmoothEdgeVertexCollision>(cc))
            {
                if (vert_edge_2_to_id.find(cc->get_hash()) == vert_edge_2_to_id.end())
                {
                    merged_collisions.collisions.push_back(cc);
                    vert_edge_2_to_id.emplace(cc->get_hash(), std::make_tuple<SmoothEdgeVertexCollision, long>(std::move(*ev), merged_collisions.collisions.size()));
                }
            }
            else if (auto ev3 = std::dynamic_pointer_cast<SmoothEdgeVertex3Collision>(cc))
            {
                if (edge_vert_3_to_id.find(cc->get_hash()) == edge_vert_3_to_id.end())
                {
                    merged_collisions.collisions.push_back(cc);
                    edge_vert_3_to_id.emplace(cc->get_hash(), std::make_tuple<SmoothEdgeVertex3Collision, long>(std::move(*ev3), merged_collisions.collisions.size()));
                }
            }
            else if (auto ff = std::dynamic_pointer_cast<SmoothFaceVertexCollision>(cc))
            {
                if (face_vert_to_id.find(cc->get_hash()) == face_vert_to_id.end())
                {
                    merged_collisions.collisions.push_back(cc);
                    face_vert_to_id.emplace(cc->get_hash(), std::make_tuple<SmoothFaceVertexCollision, long>(std::move(*ff), merged_collisions.collisions.size()));
                }
            }
            else
                throw std::runtime_error("Invalid collision type!");
        }
}

template class SmoothCollisionsBuilder<2>;
template class SmoothCollisionsBuilder<3>;

} // namespace ipc