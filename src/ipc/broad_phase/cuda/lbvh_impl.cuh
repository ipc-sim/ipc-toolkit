// Definition of the pimpl struct of ipc::cuda::LBVH. This header is CUDA-only
// and must be included from the ipc::cuda implementation files (.cu)
// exclusively.

#pragma once

#include <ipc/config.hpp>

#ifdef IPC_TOOLKIT_WITH_CUDA

#include <ipc/broad_phase/cuda/lbvh.hpp>
#include <ipc/utils/cuda/device_buffer.cuh>

#include <cstdint>
#include <mutex>
#include <vector>

namespace ipc::cuda {

struct LBVH::Impl {
    /// @brief A single BVH: the node array (root at index 0, same 32-byte
    /// ipc::LBVH::Node layout as the CPU path) plus the per-node Morton-sorted
    /// rightmost-leaf index used to skip subtrees in triangular traversal.
    struct DeviceBVH {
        DeviceBuffer<ipc::LBVH::Node> nodes;
        DeviceBuffer<int32_t> rightmost_leaves;

        /// @brief The number of leaves, derived from the node count (a BVH
        /// over n primitives has 2n - 1 nodes) so the two cannot disagree.
        int n_leaves() const
        {
            return nodes.empty() ? 0 : static_cast<int>((nodes.size() + 1) / 2);
        }

        void clear()
        {
            nodes.clear();
            rightmost_leaves.clear();
        }
    };

    DeviceBVH vertex_bvh;
    DeviceBVH edge_bvh;
    DeviceBVH face_bvh;

    /// @brief Device-resident candidate pairs (SoA) for one collision type,
    /// connectivity-filtered on the device. For the default (accept-all) vertex
    /// filter this is already the exact candidate set; otherwise it is a
    /// superset the host trims with the user filter.
    ///
    /// The buffers are grow-only pools. Their capacity is the largest count any
    /// traversal of this type on this object has needed, and it survives
    /// clear() -- which build() calls every timestep -- so a per-frame call
    /// sizes its first pass from the previous frame's count and pays the
    /// overflow-and-retry only when the count grows past every prior call.
    /// Nothing is initialized: the kernel writes exactly the slots it fills.
    /// The memory is released only by the destructor.
    struct DeviceCandidates {
        DeviceBuffer<int32_t> a;
        DeviceBuffer<int32_t> b;
        size_t count = 0; ///< The number of valid pairs in a/b.

        void clear()
        {
            count = 0;
            a.clear();
            b.clear();
        }
    };

    DeviceCandidates vv_candidates;
    DeviceCandidates ev_candidates;
    DeviceCandidates ee_candidates;
    DeviceCandidates fv_candidates;
    DeviceCandidates ef_candidates;
    DeviceCandidates ff_candidates;

    /// @brief Mesh connectivity, flat row-major (2 ids per edge, 3 per face).
    /// On the device for the traversal's shared-vertex filter; mirrored on the
    /// host for the can_*_collide filters, which run only when the user vertex
    /// filter is not accept-all.
    DeviceBuffer<int32_t> edges;
    DeviceBuffer<int32_t> faces;
    std::vector<int32_t> h_edges;
    std::vector<int32_t> h_faces;

    /// @brief The Morton-normalization domain: the union of the vertex boxes.
    /// A plain aggregate so it is trivially copyable and CUB can reduce it.
    struct Domain {
        double min[3];
        double max[3];
    };

    // -- Build scratch ------------------------------------------------------
    // Persistent and grow-only so a per-frame rebuild allocates nothing (every
    // cudaFree synchronizes the whole device). Nothing here is initialized
    // except where a kernel needs it: the construction infos are zeroed and
    // the roots set to -1 on every build.

    DeviceBuffer<double> vertices_t0; ///< Uploaded positions (column-major).
    DeviceBuffer<double> vertices_t1; ///< Uploaded positions (column-major).
    DeviceBuffer<double> vbox_min;    ///< Vertex box min corners (3 per box).
    DeviceBuffer<double> vbox_max;    ///< Vertex box max corners (3 per box).
    DeviceBuffer<double> ebox_min;    ///< Edge box min corners (3 per box).
    DeviceBuffer<double> ebox_max;    ///< Edge box max corners (3 per box).
    DeviceBuffer<double> fbox_min;    ///< Face box min corners (3 per box).
    DeviceBuffer<double> fbox_max;    ///< Face box max corners (3 per box).
    DeviceBuffer<Domain> domain;      ///< The normalization domain (1).
    /// @brief The Morton codes (keys) and box ids (values) to sort, as the two
    /// ping-pong buffers each an out-of-place radix sort needs: the codes
    /// kernel fills [0], and CUB alternates between [0] and [1] per pass (see
    /// build_tree()).
    DeviceBuffer<uint64_t> morton_codes[2];
    DeviceBuffer<int32_t> box_ids[2];
    DeviceBuffer<ipc::LBVH::ConstructionInfo<int>> construction_infos;
    DeviceBuffer<int> roots;                  ///< One root index per BVH (3).
    DeviceBuffer<unsigned char> sort_temp;    ///< CUB radix-sort storage.
    DeviceBuffer<unsigned char> reduce_temp;  ///< CUB reduce storage.
    DeviceBuffer<unsigned long long> counter; ///< Emitted-pair counter (1).

    /// @brief Serializes the detect_*() calls, which are const on the
    /// BroadPhase interface but share the candidate buffers and the counter.
    std::mutex mutex;

    /// @brief Forget the built trees, connectivity, and candidates. Every
    /// allocation is kept for the next build (see DeviceCandidates).
    void clear()
    {
        vertex_bvh.clear();
        edge_bvh.clear();
        face_bvh.clear();
        edges.clear();
        faces.clear();
        h_edges.clear();
        h_faces.clear();
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
