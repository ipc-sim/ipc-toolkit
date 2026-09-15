#pragma once

#include <ipc/collision_filter.hpp>
#include <ipc/config.hpp>

#include <cstddef>

namespace ipc::details {

// Mesh-connectivity collision filters shared by every broad phase.
//
// ipc::BroadPhase, ipc::LBVH, ipc::SweepAndTiniestQueue, and ipc::cuda::LBVH
// each store the connectivity differently -- in the AABBs' vertex_ids, in a
// dedicated host copy, or in a host mirror of the device arrays -- but they all
// apply the same rule: exclude primitive pairs that share a vertex, then accept
// the pair only if the user vertex filter passes for at least one of the
// remaining vertex pairs. These take the vertex ids directly so each broad
// phase can supply them from whatever storage it has.
//
// The rule is split in two so its halves can live where they are needed:
// share_vertex() has no dependency on the user filter and is host/device, so
// the device traversal of ipc::cuda::LBVH applies exactly the same exclusion
// the host does; any_vertex_pair_can_collide() is the user-filter half, which
// only ever runs on the host.

/// @brief Whether two primitives share a vertex id.
///
/// A vertex primitive has one id, an edge two, a face three. Fully unrolled at
/// compile time, so every id stays in a register on the device.
///
/// @tparam NA The number of vertex ids of the first primitive (1, 2, or 3).
/// @tparam NB The number of vertex ids of the second primitive (1, 2, or 3).
/// @tparam Id The vertex id type.
/// @param a The first primitive's vertex ids.
/// @param b The second primitive's vertex ids.
/// @return Whether any id of @p a equals any id of @p b.
template <int NA, int NB, typename Id>
IPC_TOOLKIT_HOST_DEVICE inline bool
share_vertex(const Id (&a)[NA], const Id (&b)[NB])
{
    static_assert(NA >= 1 && NA <= 3 && NB >= 1 && NB <= 3);
    bool shared = false;
    // The unroll pragma is nvcc's; the host compilers unroll these constant
    // trip counts on their own and would only warn about the unknown pragma.
#ifdef __CUDA_ARCH__
#pragma unroll
#endif
    for (int i = 0; i < NA; ++i) {
#ifdef __CUDA_ARCH__
#pragma unroll
#endif
        for (int j = 0; j < NB; ++j) {
            shared |= a[i] == b[j];
        }
    }
    return shared;
}

/// @brief Whether the user vertex filter passes for at least one pair of
/// vertices drawn from two primitives.
/// @tparam NA The number of vertex ids of the first primitive (1, 2, or 3).
/// @tparam NB The number of vertex ids of the second primitive (1, 2, or 3).
/// @tparam Id The vertex id type.
/// @param a The first primitive's vertex ids.
/// @param b The second primitive's vertex ids.
/// @param can_vertices_collide The user vertex filter.
/// @return Whether any vertex pair passes the filter. Always true for the
/// accept-all filter, without evaluating anything.
template <int NA, int NB, typename Id>
inline bool any_vertex_pair_can_collide(
    const Id (&a)[NA],
    const Id (&b)[NB],
    const CollisionFilter& can_vertices_collide)
{
    static_assert(NA >= 1 && NA <= 3 && NB >= 1 && NB <= 3);
    if (can_vertices_collide.accepts_all()) {
        return true;
    }
    for (int i = 0; i < NA; ++i) {
        for (int j = 0; j < NB; ++j) {
            if (can_vertices_collide(a[i], b[j])) {
                return true;
            }
        }
    }
    return false;
}

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
    const index_t e[2] = { e0i, e1i };
    const index_t v[1] = { static_cast<index_t>(vi) };
    return !share_vertex(e, v)
        && any_vertex_pair_can_collide(v, e, can_vertices_collide);
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
    const index_t ea[2] = { ea0i, ea1i };
    const index_t eb[2] = { eb0i, eb1i };
    return !share_vertex(ea, eb)
        && any_vertex_pair_can_collide(ea, eb, can_vertices_collide);
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
    const index_t f[3] = { f0i, f1i, f2i };
    const index_t v[1] = { static_cast<index_t>(vi) };
    return !share_vertex(f, v)
        && any_vertex_pair_can_collide(v, f, can_vertices_collide);
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
    const index_t e[2] = { e0i, e1i };
    const index_t f[3] = { f0i, f1i, f2i };
    return !share_vertex(e, f)
        && any_vertex_pair_can_collide(e, f, can_vertices_collide);
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
    const index_t fa[3] = { fa0i, fa1i, fa2i };
    const index_t fb[3] = { fb0i, fb1i, fb2i };
    return !share_vertex(fa, fb)
        && any_vertex_pair_can_collide(fa, fb, can_vertices_collide);
}

} // namespace ipc::details
