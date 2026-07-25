// Definition of the pimpl struct of ipc::cuda::LBVH. This header is CUDA-only
// and must be included from the ipc::cuda implementation files (.cu)
// exclusively.

#pragma once

#include <ipc/config.hpp>

#ifdef IPC_TOOLKIT_WITH_CUDA

#include <ipc/broad_phase/cuda/lbvh.hpp>

#include <thrust/device_vector.h>

#include <array>
#include <vector>

namespace ipc::cuda {

struct LBVH::Impl {
    /// @brief A single BVH: the node array (root at index 0, same 32-byte
    /// ipc::LBVH::Node layout as the CPU path) plus the per-node Morton-sorted
    /// rightmost-leaf index used to skip subtrees in triangular traversal.
    struct DeviceBVH {
        thrust::device_vector<ipc::LBVH::Node> nodes;
        thrust::device_vector<int32_t> rightmost_leaves;
        int n_leaves = 0;

        void clear()
        {
            nodes.clear();
            rightmost_leaves.clear();
            n_leaves = 0;
        }
    };

    DeviceBVH vertex_bvh;
    DeviceBVH edge_bvh;
    DeviceBVH face_bvh;

    /// @brief Device-resident candidate pairs (SoA) for one collision type,
    /// connectivity-filtered on the device. For the default (accept-all) vertex
    /// filter this is already the exact candidate set; otherwise it is a
    /// superset the host trims with the user filter.
    struct DeviceCandidates {
        thrust::device_vector<int32_t> a;
        thrust::device_vector<int32_t> b;

        /// @brief Largest candidate count ever observed for this type on this
        /// object, used to size the next traversal's output buffer so repeated
        /// calls (e.g. one per Newton iteration, or one per build() at a new
        /// timestep) don't pay the overflow-and-retry cost every time -- only
        /// the first time, or when the count grows past every prior call.
        /// Deliberately NOT reset by clear() (see below): build() calls
        /// clear() every timestep, and this hint must survive that so the
        /// learned size doesn't need re-discovering each time.
        size_t predicted_capacity = 0;

        void clear()
        {
            a.clear();
            b.clear();
            // predicted_capacity is intentionally left untouched.
        }
    };

    DeviceCandidates vv_candidates;
    DeviceCandidates ev_candidates;
    DeviceCandidates ee_candidates;
    DeviceCandidates fv_candidates;
    DeviceCandidates ef_candidates;
    DeviceCandidates ff_candidates;

    // Mesh connectivity, uploaded once and used by the device traversal's
    // shared-vertex (connectivity) filter. Flat row-major:
    // edges = 2 * n_edges, faces = 3 * n_faces.
    thrust::device_vector<index_t> edges;
    thrust::device_vector<index_t> faces;

    // Host copies of the connectivity, used by the host-side can_*_collide
    // filters applied to the device-emitted candidate pairs.
    std::vector<std::array<index_t, 2>> h_edge_vertex_ids;
    std::vector<std::array<index_t, 3>> h_face_vertex_ids;

    void clear()
    {
        vertex_bvh.clear();
        edge_bvh.clear();
        face_bvh.clear();
        edges.clear();
        faces.clear();
        h_edge_vertex_ids.clear();
        h_face_vertex_ids.clear();
        vv_candidates.clear();
        ev_candidates.clear();
        ee_candidates.clear();
        fv_candidates.clear();
        ef_candidates.clear();
        ff_candidates.clear();
    }
};

} // namespace ipc::cuda

#endif // IPC_TOOLKIT_WITH_CUDA
