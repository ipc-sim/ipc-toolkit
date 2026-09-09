#pragma once

#include <ipc/collision_filter.hpp>
#include <ipc/config.hpp>

#include <cstddef>

namespace ipc::details {

// Mesh-connectivity collision filters shared by every broad phase.
//
// ipc::BroadPhase, ipc::LBVH, and ipc::cuda::LBVH each store the connectivity
// differently -- in the AABBs' vertex_ids, in a dedicated host copy, or in a
// host mirror of the device arrays -- but they all apply the same rule: exclude
// primitive pairs that share a vertex, then accept the pair only if the user
// vertex filter passes for at least one of the remaining vertex pairs. These
// take the vertex ids directly so each broad phase can supply them from
// whatever storage it has.

/// @brief Whether an edge and a vertex can collide.
/// @param e0i The first vertex of the edge.
/// @param e1i The second vertex of the edge.
/// @param vi The vertex.
/// @param can_vertices_collide The user vertex filter.
/// @return Whether the pair should be considered for collision.
inline bool can_edge_vertex_collide(
    const index_t e0i,
    const index_t e1i,
    const size_t vi,
    const CollisionFilter& can_vertices_collide)
{
    return vi != e0i && vi != e1i
        && (can_vertices_collide(vi, e0i) || can_vertices_collide(vi, e1i));
}

/// @brief Whether two edges can collide.
/// @param ea0i The first vertex of the first edge.
/// @param ea1i The second vertex of the first edge.
/// @param eb0i The first vertex of the second edge.
/// @param eb1i The second vertex of the second edge.
/// @param can_vertices_collide The user vertex filter.
/// @return Whether the pair should be considered for collision.
inline bool can_edges_collide(
    const index_t ea0i,
    const index_t ea1i,
    const index_t eb0i,
    const index_t eb1i,
    const CollisionFilter& can_vertices_collide)
{
    const bool share_endpoint =
        ea0i == eb0i || ea0i == eb1i || ea1i == eb0i || ea1i == eb1i;

    return !share_endpoint
        && (can_vertices_collide(ea0i, eb0i) || can_vertices_collide(ea0i, eb1i)
            || can_vertices_collide(ea1i, eb0i)
            || can_vertices_collide(ea1i, eb1i));
}

/// @brief Whether a face and a vertex can collide.
/// @param f0i The first vertex of the face.
/// @param f1i The second vertex of the face.
/// @param f2i The third vertex of the face.
/// @param vi The vertex.
/// @param can_vertices_collide The user vertex filter.
/// @return Whether the pair should be considered for collision.
inline bool can_face_vertex_collide(
    const index_t f0i,
    const index_t f1i,
    const index_t f2i,
    const size_t vi,
    const CollisionFilter& can_vertices_collide)
{
    return vi != f0i && vi != f1i && vi != f2i
        && (can_vertices_collide(vi, f0i) || can_vertices_collide(vi, f1i)
            || can_vertices_collide(vi, f2i));
}

/// @brief Whether an edge and a face can intersect.
/// @param e0i The first vertex of the edge.
/// @param e1i The second vertex of the edge.
/// @param f0i The first vertex of the face.
/// @param f1i The second vertex of the face.
/// @param f2i The third vertex of the face.
/// @param can_vertices_collide The user vertex filter.
/// @return Whether the pair should be considered for collision.
inline bool can_edge_face_collide(
    const index_t e0i,
    const index_t e1i,
    const index_t f0i,
    const index_t f1i,
    const index_t f2i,
    const CollisionFilter& can_vertices_collide)
{
    const bool share_endpoint = e0i == f0i || e0i == f1i || e0i == f2i
        || e1i == f0i || e1i == f1i || e1i == f2i;

    return !share_endpoint
        && (can_vertices_collide(e0i, f0i) || can_vertices_collide(e0i, f1i)
            || can_vertices_collide(e0i, f2i) || can_vertices_collide(e1i, f0i)
            || can_vertices_collide(e1i, f1i)
            || can_vertices_collide(e1i, f2i));
}

/// @brief Whether two faces can collide.
/// @param fa0i The first vertex of the first face.
/// @param fa1i The second vertex of the first face.
/// @param fa2i The third vertex of the first face.
/// @param fb0i The first vertex of the second face.
/// @param fb1i The second vertex of the second face.
/// @param fb2i The third vertex of the second face.
/// @param can_vertices_collide The user vertex filter.
/// @return Whether the pair should be considered for collision.
inline bool can_faces_collide(
    const index_t fa0i,
    const index_t fa1i,
    const index_t fa2i,
    const index_t fb0i,
    const index_t fb1i,
    const index_t fb2i,
    const CollisionFilter& can_vertices_collide)
{
    const bool share_endpoint = fa0i == fb0i || fa0i == fb1i || fa0i == fb2i
        || fa1i == fb0i || fa1i == fb1i || fa1i == fb2i || fa2i == fb0i
        || fa2i == fb1i || fa2i == fb2i;

    return !share_endpoint
        && (can_vertices_collide(fa0i, fb0i) //
            || can_vertices_collide(fa0i, fb1i)
            || can_vertices_collide(fa0i, fb2i)
            || can_vertices_collide(fa1i, fb0i)
            || can_vertices_collide(fa1i, fb1i)
            || can_vertices_collide(fa1i, fb2i)
            || can_vertices_collide(fa2i, fb0i)
            || can_vertices_collide(fa2i, fb1i)
            || can_vertices_collide(fa2i, fb2i));
}

} // namespace ipc::details
