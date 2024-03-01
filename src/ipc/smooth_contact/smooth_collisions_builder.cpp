#include "smooth_collisions_builder.hpp"
#include <tbb/enumerable_thread_specific.h>
#include <ipc/distance/distance_type.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/point_triangle.hpp>
#include <iostream>

namespace ipc {

namespace {
    template <int dim, typename TCollision>
    void add_collision(
        const std::shared_ptr<TCollision>& pair,
        unordered_map<std::pair<long, long>, std::shared_ptr<TCollision>>&
            cc_to_id_,
        std::vector<
            std::shared_ptr<typename SmoothCollisions<dim>::value_type>>&
            collisions_)
    {
        if (pair->is_active()
            && cc_to_id_.find(pair->get_hash()) == cc_to_id_.end()) {
            // New collision, so add it to the end of collisions
            cc_to_id_.emplace(
                pair->get_hash(),
                pair);
            collisions_.push_back(pair);
        }
    }
} // namespace

void SmoothCollisionsBuilder<2>::add_edge_vertex_collisions(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const std::vector<EdgeVertexCandidate>& candidates,
    const ParameterType& param,
    const std::function<double(const long&)>& vert_dhat,
    const std::function<double(const long&)>& edge_dhat,
    const size_t start_i,
    const size_t end_i)
{
    for (size_t i = start_i; i < end_i; i++) {
        const auto& [ei, vi] = candidates[i];

        add_collision<2, SmoothCollisionTemplate<max_vert_2d, Edge2, Point2>>(
            std::make_shared<
                SmoothCollisionTemplate<max_vert_2d, Edge2, Point2>>(
                ei, vi, PointEdgeDistanceType::AUTO, mesh, param,
                std::min(edge_dhat(ei), vert_dhat(vi)), vertices),
            vert_edge_2_to_id, collisions);

        for (int j : { 0, 1 }) {
            add_collision<
                2, SmoothCollisionTemplate<max_vert_2d, Point2, Point2>>(
                std::make_shared<
                    SmoothCollisionTemplate<max_vert_2d, Point2, Point2>>(
                    mesh.edges()(ei, j), vi, PointPointDistanceType::AUTO, mesh,
                    param,
                    std::min(vert_dhat(mesh.edges()(ei, j)), vert_dhat(vi)),
                    vertices),
                vert_vert_2_to_id, collisions);
        }
    }
}

// ============================================================================

void SmoothCollisionsBuilder<3>::add_edge_edge_collisions(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const std::vector<EdgeEdgeCandidate>& candidates,
    const ParameterType& param,
    const std::function<double(const long&)>& vert_dhat,
    const std::function<double(const long&)>& edge_dhat,
    const size_t start_i,
    const size_t end_i)
{
    for (size_t i = start_i; i < end_i; i++) {
        const auto& [eai, ebi] = candidates[i];

        const auto [ea0, ea1, eb0, eb1] =
            candidates[i].vertices(vertices, mesh.edges(), mesh.faces());

        const EdgeEdgeDistanceType actual_dtype =
            edge_edge_distance_type(ea0, ea1, eb0, eb1);

        const double distance =
            sqrt(edge_edge_distance(ea0, ea1, eb0, eb1, actual_dtype));

        if (actual_dtype != EdgeEdgeDistanceType::EA_EB
            || distance >= param.dhat)
            continue;

        add_collision<3, SmoothCollisionTemplate<max_vert_3d, Edge3, Edge3>>(
            std::make_shared<
                SmoothCollisionTemplate<max_vert_3d, Edge3, Edge3>>(
                std::min(eai, ebi), std::max(eai, ebi), actual_dtype, mesh,
                param, std::min(edge_dhat(eai), edge_dhat(ebi)), vertices),
            edge_edge_3_to_id, collisions);
    }
}

void SmoothCollisionsBuilder<3>::add_face_vertex_collisions(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const std::vector<FaceVertexCandidate>& candidates,
    const ParameterType& param,
    const std::function<double(const long&)>& vert_dhat,
    const std::function<double(const long&)>& edge_dhat,
    const std::function<double(const long&)>& face_dhat,
    const size_t start_i,
    const size_t end_i)
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
            add_collision<
                3, SmoothCollisionTemplate<max_vert_3d, Face, Point3>>(
                std::make_shared<
                    SmoothCollisionTemplate<max_vert_3d, Face, Point3>>(
                    fi, vi, pt_dtype, mesh, param,
                    std::min(face_dhat(fi), vert_dhat(vi)), vertices),
                face_vert_to_id, collisions);

        for (int lv = 0; lv < 3; lv++) {
            const auto& vj = mesh.faces()(fi, lv);
            const double dhat = std::min(vert_dhat(vi), vert_dhat(vj));
            if ((vertices.row(vi) - vertices.row(vj)).norm() >= dhat)
                continue;
            add_collision<
                3, SmoothCollisionTemplate<max_vert_3d, Point3, Point3>>(
                std::make_shared<
                    SmoothCollisionTemplate<max_vert_3d, Point3, Point3>>(
                    std::min<long>(vi, vj), std::max<long>(vi, vj),
                    PointPointDistanceType::AUTO, mesh, param, dhat, vertices),
                vert_vert_3_to_id, collisions);
        }

        for (int le = 0; le < 3; le++) {
            const auto& eid = mesh.faces_to_edges()(fi, le);
            const double dhat = std::min(edge_dhat(eid), vert_dhat(vi));

            const PointEdgeDistanceType pe_dtype = point_edge_distance_type(
                vertices.row(vi), vertices.row(mesh.edges()(eid, 0)),
                vertices.row(mesh.edges()(eid, 1)));
            const double distance_sqr = point_edge_distance(
                vertices.row(vi), vertices.row(mesh.edges()(eid, 0)),
                vertices.row(mesh.edges()(eid, 1)), pe_dtype);

            if (pe_dtype != PointEdgeDistanceType::P_E
                || sqrt(distance_sqr) >= dhat)
                continue;

            add_collision<
                3, SmoothCollisionTemplate<max_vert_3d, Edge3, Point3>>(
                std::make_shared<
                    SmoothCollisionTemplate<max_vert_3d, Edge3, Point3>>(
                    eid, vi, pe_dtype, mesh, param, dhat, vertices),
                edge_vert_3_to_id, collisions);
        }
    }
}

void SmoothCollisionsBuilder<3>::merge(
    const utils::ParallelCacheType<SmoothCollisionsBuilder<3>>& local_storage,
    SmoothCollisions<3>& merged_collisions)
{
    unordered_map<
        std::pair<long, long>,
        std::shared_ptr<SmoothCollisionTemplate<max_vert_3d, Face, Point3>>>
        face_vert_to_id;
    unordered_map<
        std::pair<long, long>,
        std::shared_ptr<SmoothCollisionTemplate<max_vert_3d, Point3, Point3>>>
        vert_vert_3_to_id;
    unordered_map<
        std::pair<long, long>,
        std::shared_ptr<SmoothCollisionTemplate<max_vert_3d, Edge3, Point3>>>
        edge_vert_3_to_id;
    unordered_map<
        std::pair<long, long>,
        std::shared_ptr<SmoothCollisionTemplate<max_vert_3d, Edge3, Edge3>>>
        edge_edge_3_to_id;

    // size up the hash items
    size_t total = 0;
    for (const auto& storage : local_storage)
        total += storage.collisions.size();

    merged_collisions.collisions.reserve(total);

    // merge
    for (auto& builder : local_storage)
    {
        face_vert_to_id.insert(builder.face_vert_to_id.begin(), builder.face_vert_to_id.end());
        vert_vert_3_to_id.insert(builder.vert_vert_3_to_id.begin(), builder.vert_vert_3_to_id.end());
        edge_vert_3_to_id.insert(builder.edge_vert_3_to_id.begin(), builder.edge_vert_3_to_id.end());
        edge_edge_3_to_id.insert(builder.edge_edge_3_to_id.begin(), builder.edge_edge_3_to_id.end());
    }
    int edge_vert_count = edge_vert_3_to_id.size();
    int vert_vert_count = vert_vert_3_to_id.size();
    int face_vert_count = face_vert_to_id.size();
    int edge_edge_count = edge_edge_3_to_id.size();

    for (const auto& [key, val] : face_vert_to_id)
        merged_collisions.collisions.push_back(val);
    for (const auto& [key, val] : vert_vert_3_to_id)
        merged_collisions.collisions.push_back(val);
    for (const auto& [key, val] : edge_vert_3_to_id)
        merged_collisions.collisions.push_back(val);
    for (const auto& [key, val] : edge_edge_3_to_id)
        merged_collisions.collisions.push_back(val);

    logger().trace(
        "edge-vert pairs {}, vert-vert pairs {}", edge_vert_count,
        vert_vert_count);
    logger().trace(
        "face-vert pairs {}, edge-edge pairs {}", face_vert_count,
        edge_edge_count);
}

void SmoothCollisionsBuilder<2>::merge(
    const utils::ParallelCacheType<SmoothCollisionsBuilder<2>>& local_storage,
    SmoothCollisions<2>& merged_collisions)
{
    unordered_map<
        std::pair<long, long>,
        std::shared_ptr<SmoothCollisionTemplate<max_vert_2d, Point2, Point2>>>
        vert_vert_2_to_id;
    unordered_map<
        std::pair<long, long>,
        std::shared_ptr<SmoothCollisionTemplate<max_vert_2d, Edge2, Point2>>>
        vert_edge_2_to_id;

    // size up the hash items
    size_t total = 0;
    for (const auto& storage : local_storage)
        total += storage.collisions.size();

    merged_collisions.collisions.reserve(total);

    // merge
    for (auto& builder : local_storage)
    {
        vert_vert_2_to_id.insert(builder.vert_vert_2_to_id.begin(), builder.vert_vert_2_to_id.end());
        vert_edge_2_to_id.insert(builder.vert_edge_2_to_id.begin(), builder.vert_edge_2_to_id.end());
    }
    int edge_vert_count = vert_edge_2_to_id.size();
    int vert_vert_count = vert_vert_2_to_id.size();

    for (const auto& [key, val] : vert_vert_2_to_id)
        merged_collisions.collisions.push_back(val);
    for (const auto& [key, val] : vert_edge_2_to_id)
        merged_collisions.collisions.push_back(val);

    logger().trace(
        "edge-vert pairs {}, vert-vert pairs {}", edge_vert_count,
        vert_vert_count);
}

template class SmoothCollisionsBuilder<2>;
template class SmoothCollisionsBuilder<3>;

} // namespace ipc