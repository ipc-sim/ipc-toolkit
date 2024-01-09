#include "smooth_collisions_builder.hpp"

#include <ipc/distance/distance_type.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/edge_edge_mollifier.hpp>
#include <ipc/distance/point_triangle.hpp>

#include <iostream>

namespace ipc {

void SmoothCollisionsBuilder::add_edge_vertex_collisions(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const std::vector<EdgeVertexCandidate>& candidates,
    const std::function<bool(double)>& is_active,
    const size_t start_i,
    const size_t end_i,
    const double dhat,
    const bool use_adaptive_eps)
{
    std::cout << "adding edge vertex collisions" << std::endl;
    for (size_t i = start_i; i < end_i; i++) {
        const auto& [ei, vi] = candidates[i];
        const auto [v, e0, e1, _] =
            candidates[i].vertices(vertices, mesh.edges(), mesh.faces());
        const PointEdgeDistanceType dtype = point_edge_distance_type(v, e0, e1);
        const double distance_sqr = point_edge_distance(v, e0, e1, dtype);

        if (!is_active(distance_sqr))
            continue;

        const double weight = mesh.vertex_area(vi) / 2;
        Eigen::SparseVector<double> weight_gradient;

        add_edge_vertex_collision(
            mesh, candidates[i], dtype, weight, weight_gradient, dhat, use_adaptive_eps);
    }
}

void SmoothCollisionsBuilder::add_edge_vertex_collision(
    const CollisionMesh& mesh,
    const EdgeVertexCandidate& candidate,
    const PointEdgeDistanceType dtype,
    const double weight,
    const Eigen::SparseVector<double>& weight_gradient,
    const double dhat,
    const bool use_adaptive_eps)
{
    const auto& [ei, vi] = candidate;
    if (use_adaptive_eps) {
        std::cout << "new dhat" << std::endl;
        add_edge_vertex_collision(ei, vi, weight, weight_gradient, std::min(dhat, mesh.min_distance_in_rest_config(vi)));
    } else {
        std::cout << "old dhat" << std::endl;
        add_edge_vertex_collision(ei, vi, weight, weight_gradient, dhat*dhat);
 
    }
}

// ============================================================================

void SmoothCollisionsBuilder::add_edge_vertex_collision(
    const SmoothEdgeVertexCollision& ev_collision,
    unordered_map<SmoothEdgeVertexCollision, long>& ev_to_id,
    std::vector<SmoothEdgeVertexCollision>& ev_collisions)
{
    auto found_item = ev_to_id.find(ev_collision);
    if (found_item != ev_to_id.end()) {
        // collision already exists, so increase weight
        ev_collisions[found_item->second].weight += ev_collision.weight;
        ev_collisions[found_item->second].weight_gradient +=
            ev_collision.weight_gradient;
    } else {
        // New collision, so add it to the end of ev_collisions
        ev_to_id.emplace(ev_collision, ev_collisions.size());
        ev_collisions.push_back(ev_collision);
    }
}

void SmoothCollisionsBuilder::add_edge_edge_collisions(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const std::vector<EdgeEdgeCandidate>& candidates,
    const std::function<bool(double)>& is_active,
    const size_t start_i,
    const size_t end_i,
    const SurfaceQuadratureType quad_type,
    const double dhat,
    const bool use_adaptive_eps)
{
    std::cout << "adding edge edge collisions" << std::endl;
    for (size_t i = start_i; i < end_i; i++) {
        const auto& [eai, ebi] = candidates[i];

        const auto [ea0i, ea1i, eb0i, eb1i] =
            candidates[i].vertex_ids(mesh.edges(), mesh.faces());

        const auto [ea0, ea1, eb0, eb1] =
            candidates[i].vertices(vertices, mesh.edges(), mesh.faces());

        const double distance_sqr =
            mesh.dim() == 2 ? edge_edge_distance_2d(ea0, ea1, eb0, eb1) : edge_edge_distance(ea0, ea1, eb0, eb1, EdgeEdgeDistanceType::AUTO);

        if (!is_active(distance_sqr))
            continue;

        const double eps_x = mesh.dim() == 2 ? 0 : edge_edge_mollifier_threshold(
            mesh.rest_positions().row(ea0i), mesh.rest_positions().row(ea1i),
            mesh.rest_positions().row(eb0i), mesh.rest_positions().row(eb1i));

        const double weight = 1;
        Eigen::SparseVector<double> weight_gradient;

        if (mesh.dim() == 2 && quad_type == SurfaceQuadratureType::SinglePoint)
        {
            const double eps0 = use_adaptive_eps ? std::min(mesh.min_distance_in_rest_config(eb0i), dhat*dhat) : dhat*dhat;
            const double eps1 = use_adaptive_eps ? std::min(mesh.min_distance_in_rest_config(eb1i), dhat*dhat) : dhat;
            const double eps2 = use_adaptive_eps ? std::min(mesh.min_distance_in_rest_config(ea0i), dhat*dhat) : dhat*dhat;
            const double eps3 = use_adaptive_eps ? std::min(mesh.min_distance_in_rest_config(ea1i), dhat*dhat) : dhat*dhat;

            add_edge_vertex_collision(eai, eb0i, mesh.vertex_area(eb0i) / 2, weight_gradient, 
                eps0);
            add_edge_vertex_collision(eai, eb1i, mesh.vertex_area(eb1i) / 2, weight_gradient, 
                eps1);

            add_edge_vertex_collision(ebi, ea0i, mesh.vertex_area(ea0i) / 2, weight_gradient, 
                eps2);
            add_edge_vertex_collision(ebi, ea1i, mesh.vertex_area(ea1i) / 2, weight_gradient, 
                eps3);
        }
        else
            add_edge_edge_collision(SmoothEdgeEdgeCollision(
                                eai, ebi, eps_x, weight, weight_gradient, EdgeEdgeDistanceType::AUTO),
                                ee_to_id, ee_collisions);
    }
}

void SmoothCollisionsBuilder::add_face_vertex_collisions(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const std::vector<FaceVertexCandidate>& candidates,
    const std::function<bool(double)>& is_active,
    const size_t start_i,
    const size_t end_i)
{
    for (size_t i = start_i; i < end_i; i++) {
        const auto& [fi, vi] = candidates[i];
        const long f0i = mesh.faces()(fi, 0), f1i = mesh.faces()(fi, 1),
                   f2i = mesh.faces()(fi, 2);

        const auto [v, f0, f1, f2] =
            candidates[i].vertices(vertices, mesh.edges(), mesh.faces());

        // Compute distance type
        const PointTriangleDistanceType dtype =
            point_triangle_distance_type(v, f0, f1, f2);
        const double distance_sqr =
            point_triangle_distance(v, f0, f1, f2, dtype);

        if (!is_active(distance_sqr))
            continue;

        const double weight = mesh.vertex_area(vi) / 3;
        Eigen::SparseVector<double> weight_gradient;
        
        fv_collisions.emplace_back(fi, vi, weight, weight_gradient);
    }
}

void SmoothCollisionsBuilder::add_edge_edge_collision(
    const SmoothEdgeEdgeCollision& ee_collision,
    unordered_map<SmoothEdgeEdgeCollision, long>& ee_to_id_,
    std::vector<SmoothEdgeEdgeCollision>& ee_collisions_)
{
    auto found_item = ee_to_id_.find(ee_collision);
    if (found_item != ee_to_id_.end()) {
        // edge-edge collision shouldn't be counted multiple times
        // assert(ee_collision == ee_collisions_[found_item->second]);
        // ee_collisions_[found_item->second].weight += ee_collision.weight;
        // ee_collisions_[found_item->second].weight_gradient +=
        //     ee_collision.weight_gradient;
    } else {
        // New collision, so add it to the end of ee_collisions
        ee_to_id_.emplace(ee_collision, ee_collisions_.size());
        ee_collisions_.push_back(ee_collision);
    }
}

void SmoothCollisionsBuilder::merge(
    const tbb::enumerable_thread_specific<SmoothCollisionsBuilder>& local_storage,
    SmoothCollisions& merged_collisions)
{
    unordered_map<SmoothEdgeVertexCollision, long> ev_to_id;
    unordered_map<SmoothEdgeEdgeCollision, long> ee_to_id;
    auto& ev_collisions = merged_collisions.ev_collisions;
    auto& ee_collisions = merged_collisions.ee_collisions;
    auto& fv_collisions = merged_collisions.fv_collisions;

    // size up the hash items
    size_t n_vv = 0, n_ev = 0, n_ee = 0, n_fv = 0;
    for (const auto& storage : local_storage) {
        // This is an conservative estimate
        n_ev += storage.ev_collisions.size();
        n_ee += storage.ee_collisions.size();
        n_fv += storage.fv_collisions.size();
    }
    ev_collisions.reserve(n_ev);
    ee_collisions.reserve(n_ee);
    fv_collisions.reserve(n_fv);

    // merge
    for (const auto& builder : local_storage) {
        if (ev_collisions.empty()) {
            ev_to_id = builder.ev_to_id;
            ev_collisions = builder.ev_collisions;
        } else {
            for (const auto& ev : builder.ev_collisions) {
                add_edge_vertex_collision(ev, ev_to_id, ev_collisions);
            }
        }

        if (ee_collisions.empty()) {
            ee_to_id = builder.ee_to_id;
            ee_collisions = builder.ee_collisions;
        } else {
            for (const auto& ee : builder.ee_collisions) {
                add_edge_edge_collision(ee, ee_to_id, ee_collisions);
            }
        }

        fv_collisions.insert(
            fv_collisions.end(), builder.fv_collisions.begin(),
            builder.fv_collisions.end());
    }

    // If positive and negative vertex-vertex collisions cancel out, remove
    // them. This can happen when edge-vertex collisions reduce to
    // vertex-vertex collisions. This will avoid unnecessary computation.
    // Same for edge-vertex collisions.
    ev_collisions.erase(
        std::remove_if(
            ev_collisions.begin(), ev_collisions.end(),
            [&](const SmoothEdgeVertexCollision& ev) { return ev.weight == 0; }),
        ev_collisions.end());
    // Same for edge-edge collisions.
    ee_collisions.erase(
        std::remove_if(
            ee_collisions.begin(), ee_collisions.end(),
            [&](const SmoothEdgeEdgeCollision& ee) { return ee.weight == 0; }),
        ee_collisions.end());
}

} // namespace ipc