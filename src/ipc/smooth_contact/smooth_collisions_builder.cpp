#include "smooth_collisions_builder.hpp"

#include <ipc/distance/distance_type.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/point_triangle.hpp>
#include <iostream>
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
        for (size_t i = start_i; i < end_i; i++) {
            const auto& [ei, vi] = candidates[i];

            add_collision<SmoothEdgeVertexCollision>(std::make_shared<SmoothEdgeVertexCollision>(ei, vi, mesh, param, std::min(edge_dhat(ei), vert_dhat(vi)), vertices), vert_edge_2_to_id, collisions);

            for (int j : {0, 1})
            {
                add_collision<SmoothVertexVertexCollision>(std::make_shared<SmoothVertexVertexCollision>(mesh.edges()(ei, j), vi, mesh, param, std::min(vert_dhat(mesh.edges()(ei, j)), vert_dhat(vi)), vertices), vert_vert_2_to_id, collisions);
            }
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
    const std::function<double(const long &)> &vert_dhat,
    const std::function<double(const long &)> &edge_dhat,
    const size_t start_i,
    const size_t end_i)
{
    if constexpr (dim == 3)
    {
        for (size_t i = start_i; i < end_i; i++) {
            const auto& [eai, ebi] = candidates[i];

            const auto [ea0, ea1, eb0, eb1] =
                candidates[i].vertices(vertices, mesh.edges(), mesh.faces());

            const EdgeEdgeDistanceType actual_dtype =
                edge_edge_distance_type(ea0, ea1, eb0, eb1);

            const double distance =
                sqrt(edge_edge_distance(ea0, ea1, eb0, eb1, actual_dtype));
            
            if (distance < 1e-12)
                logger().warn("edge {} {} dist {}", eai, ebi, distance);
            
            if (actual_dtype != EdgeEdgeDistanceType::EA_EB || distance >= param.dhat)
                continue;

            add_collision<SmoothCollisionTemplate<max_vert_3d, Edge3 , Edge3 >>(std::make_shared<SmoothCollisionTemplate<max_vert_3d, Edge3 , Edge3 >>(std::min(eai, ebi), std::max(eai, ebi), actual_dtype, mesh, param, std::min(edge_dhat(eai), edge_dhat(ebi)), vertices), edge_edge_3_to_id, collisions);
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
        for (size_t i = start_i; i < end_i; i++) {
            const auto& [fi, vi] = candidates[i];

            Eigen::Vector3d v = vertices.row(vi);
            Eigen::Vector3d f0 = vertices.row(mesh.faces()(fi, 0));
            Eigen::Vector3d f1 = vertices.row(mesh.faces()(fi, 1));
            Eigen::Vector3d f2 = vertices.row(mesh.faces()(fi, 2));
            const PointTriangleDistanceType pt_dtype =
                point_triangle_distance_type(v, f0, f1, f2);
            const double distance =
                sqrt(point_triangle_distance(v, f0, f1, f2, pt_dtype));
            
            if (distance >= vert_dhat(vi))
                continue;

            if (pt_dtype == PointTriangleDistanceType::P_T)
                add_collision<SmoothCollisionTemplate<max_vert_3d, Face  , Point3>>(std::make_shared<SmoothCollisionTemplate<max_vert_3d, Face  , Point3>>(fi, vi, pt_dtype, mesh, param, std::min(face_dhat(fi), vert_dhat(vi)), vertices), face_vert_to_id, collisions);
            
            for (int lv = 0; lv < 3; lv++)
            {
                const auto &vj = mesh.faces()(fi, lv);
                const double dhat = std::min(vert_dhat(vi), vert_dhat(vj));
                if ((vertices.row(vi) - vertices.row(vj)).norm() >= dhat)
                    continue;
                add_collision<SmoothCollisionTemplate<max_vert_3d, Point3, Point3>>(std::make_shared<SmoothCollisionTemplate<max_vert_3d, Point3, Point3>>(std::min<long>(vi, vj), std::max<long>(vi, vj), PointPointDistanceType::AUTO, mesh, param, dhat, vertices), vert_vert_3_to_id, collisions);
            }
            
            for (int le = 0; le < 3; le++)
            {
                const auto &eid = mesh.faces_to_edges()(fi, le);
                const double dhat = std::min(edge_dhat(eid), vert_dhat(vi));

                const PointEdgeDistanceType pe_dtype = point_edge_distance_type(vertices.row(vi), vertices.row(mesh.edges()(eid, 0)), vertices.row(mesh.edges()(eid, 1)));
                const double distance_sqr = point_edge_distance(vertices.row(vi), vertices.row(mesh.edges()(eid, 0)), vertices.row(mesh.edges()(eid, 1)), pe_dtype);

                if (pe_dtype != PointEdgeDistanceType::P_E || sqrt(distance_sqr) >= dhat)
                    continue;

                add_collision<SmoothCollisionTemplate<max_vert_3d, Edge3 , Point3>>(std::make_shared<SmoothCollisionTemplate<max_vert_3d, Edge3 , Point3>>(eid, vi, pe_dtype, mesh, param, dhat, vertices), edge_vert_3_to_id, collisions);
            }
        }
    }
}

template <int dim> template <typename TCollision>
void SmoothCollisionsBuilder<dim>::add_collision(
    const std::shared_ptr<TCollision>& pair,
    unordered_map<std::pair<long, long>, std::tuple<TCollision, long> >& cc_to_id_,
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
    unordered_map<std::pair<long, long>, std::tuple<SmoothVertexVertexCollision, long> > vert_vert_2_to_id;
    unordered_map<std::pair<long, long>, std::tuple<SmoothEdgeVertexCollision, long> > vert_edge_2_to_id;
    
    unordered_map<std::pair<long, long>, std::tuple<SmoothCollisionTemplate<max_vert_3d, Face  , Point3>, long> > face_vert_to_id;
    unordered_map<std::pair<long, long>, std::tuple<SmoothCollisionTemplate<max_vert_3d, Point3, Point3>, long> > vert_vert_3_to_id;
    unordered_map<std::pair<long, long>, std::tuple<SmoothCollisionTemplate<max_vert_3d, Edge3 , Point3>, long> > edge_vert_3_to_id;
    unordered_map<std::pair<long, long>, std::tuple<SmoothCollisionTemplate<max_vert_3d, Edge3 , Edge3 >, long> > edge_edge_3_to_id;

    // size up the hash items
    size_t total = 0;
    for (const auto& storage : local_storage)
        total += storage.collisions.size();
    
    merged_collisions.collisions.reserve(total);

    int edge_vert_count = 0;
    int vert_vert_count = 0;
    int face_vert_count = 0;
    int edge_edge_count = 0;

    // merge
    for (auto& builder : local_storage)
        for (auto& cc : builder.collisions)
        {
            if (auto ee3 = std::dynamic_pointer_cast<SmoothCollisionTemplate<max_vert_3d, Edge3 , Edge3 >>(cc))
            {
                if (edge_edge_3_to_id.find(cc->get_hash()) == edge_edge_3_to_id.end())
                {
                    merged_collisions.collisions.push_back(cc);
                    edge_edge_3_to_id.emplace(cc->get_hash(), std::make_tuple<SmoothCollisionTemplate<max_vert_3d, Edge3 , Edge3 >, long>(std::move(*ee3), merged_collisions.collisions.size()));
                    edge_edge_count++;
                }
            }
            else if (auto vv2 = std::dynamic_pointer_cast<SmoothVertexVertexCollision>(cc))
            {
                if (vert_vert_2_to_id.find(cc->get_hash()) == vert_vert_2_to_id.end())
                {
                    merged_collisions.collisions.push_back(cc);
                    vert_vert_2_to_id.emplace(cc->get_hash(), std::make_tuple<SmoothVertexVertexCollision, long>(std::move(*vv2), merged_collisions.collisions.size()));
                    vert_vert_count++;
                }
            }
            else if (auto vv3 = std::dynamic_pointer_cast<SmoothCollisionTemplate<max_vert_3d, Point3, Point3>>(cc))
            {
                if (vert_vert_3_to_id.find(cc->get_hash()) == vert_vert_3_to_id.end())
                {
                    merged_collisions.collisions.push_back(cc);
                    vert_vert_3_to_id.emplace(cc->get_hash(), std::make_tuple<SmoothCollisionTemplate<max_vert_3d, Point3, Point3>, long>(std::move(*vv3), merged_collisions.collisions.size()));
                    vert_vert_count++;
                }
            }
            else if (auto ev = std::dynamic_pointer_cast<SmoothEdgeVertexCollision>(cc))
            {
                if (vert_edge_2_to_id.find(cc->get_hash()) == vert_edge_2_to_id.end())
                {
                    merged_collisions.collisions.push_back(cc);
                    vert_edge_2_to_id.emplace(cc->get_hash(), std::make_tuple<SmoothEdgeVertexCollision, long>(std::move(*ev), merged_collisions.collisions.size()));
                    edge_vert_count++;
                }
            }
            else if (auto ev3 = std::dynamic_pointer_cast<SmoothCollisionTemplate<max_vert_3d, Edge3 , Point3>>(cc))
            {
                if (edge_vert_3_to_id.find(cc->get_hash()) == edge_vert_3_to_id.end())
                {
                    merged_collisions.collisions.push_back(cc);
                    edge_vert_3_to_id.emplace(cc->get_hash(), std::make_tuple<SmoothCollisionTemplate<max_vert_3d, Edge3 , Point3>, long>(std::move(*ev3), merged_collisions.collisions.size()));
                    edge_vert_count++;
                }
            }
            else if (auto ff = std::dynamic_pointer_cast<SmoothCollisionTemplate<max_vert_3d, Face  , Point3>>(cc))
            {
                if (face_vert_to_id.find(cc->get_hash()) == face_vert_to_id.end())
                {
                    merged_collisions.collisions.push_back(cc);
                    face_vert_to_id.emplace(cc->get_hash(), std::make_tuple<SmoothCollisionTemplate<max_vert_3d, Face  , Point3>, long>(std::move(*ff), merged_collisions.collisions.size()));
                    face_vert_count++;
                }
            }
            else
                throw std::runtime_error("Invalid collision type!");
        }

    logger().trace("edge-vert pairs {}, vert-vert pairs {}", edge_vert_count, vert_vert_count);
    if (face_vert_count || edge_edge_count)
        logger().trace("face-vert pairs {}, edge-edge pairs {}", face_vert_count, edge_edge_count);
}

template class SmoothCollisionsBuilder<2>;
template class SmoothCollisionsBuilder<3>;

} // namespace ipc