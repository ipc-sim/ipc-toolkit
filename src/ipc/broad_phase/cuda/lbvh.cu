#include "lbvh.hpp"

#ifdef IPC_TOOLKIT_WITH_CUDA

#include <ipc/broad_phase/cuda/lbvh_impl.cuh>
#include <ipc/broad_phase/details/connectivity_filters.hpp>
#include <ipc/broad_phase/details/lbvh_build.hpp>
#include <ipc/broad_phase/details/lbvh_traverse.hpp>
#include <ipc/math/morton.hpp>
#include <ipc/utils/cuda/device_utils.cuh>
#include <ipc/utils/logger.hpp>

#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>

#include <algorithm>
#include <cassert>
#include <cub/device/device_radix_sort.cuh>
#include <cub/device/device_reduce.cuh>
#include <limits>
#include <mutex>
#include <string>
#include <type_traits>
#include <vector>

namespace ipc::cuda {

namespace {

    // ipc::LBVH::Node is shared through device memory and copied back to the
    // host bytewise, so both sides must agree on its layout. ipc::LBVH asserts
    // the size in a host-only constructor no device path instantiates; this
    // one is evaluated by nvcc's front end, which is what lays the type out for
    // the device.
    static_assert(
        sizeof(ipc::LBVH::Node) == 32,
        "ipc::LBVH::Node must be 32 bytes to share it with the host");
    static_assert(
        alignof(ipc::LBVH::Node) == 32,
        "ipc::LBVH::Node must be 32-byte aligned to share it with the host");

    /// @brief Per-internal-node scratch used by the bottom-up build. The
    /// counter is a plain int rather than the host build's std::atomic<int>:
    /// atomicAdd() needs an int*, and supplies the same atomicity. The layout
    /// is otherwise identical to the host's.
    using DeviceConstructionInfo = ipc::LBVH::ConstructionInfo<int>;

    using Domain = LBVH::Impl::Domain;
    static_assert(
        std::is_trivially_copyable_v<Domain>,
        "Domain is copied into kernel parameter space and reduced by CUB");

    struct DomainReduce {
        __host__ __device__ Domain
        operator()(const Domain& a, const Domain& b) const
        {
            Domain r;
#pragma unroll
            for (int k = 0; k < 3; ++k) {
                r.min[k] = fmin(a.min[k], b.min[k]);
                r.max[k] = fmax(a.max[k], b.max[k]);
            }
            return r;
        }
    };

    struct MakeDomain {
        const double* box_min;
        const double* box_max;
        __host__ __device__ Domain operator()(const int i) const
        {
            Domain d;
#pragma unroll
            for (int k = 0; k < 3; ++k) {
                d.min[k] = box_min[3 * i + k];
                d.max[k] = box_max[3 * i + k];
            }
            return d;
        }
    };

    // -- Box building -------------------------------------------------------

    /// @brief One vertex box from the vertex's positions at t0 and t1; pass
    /// the same array twice for a static box.
    ///
    /// Mirrors AABB::from_point(p_t0, p_t1, r), i.e. the union of the two
    /// inflated points, bit-for-bit: the per-coordinate inflation is
    /// AABB::conservative_{lower,upper}_bound() -- the very function the host
    /// calls -- and taking the union before inflating equals inflating before
    /// the union because nextafter is monotone:
    /// min(nextafter(a - r), nextafter(b - r)) == nextafter(min(a, b) - r).
    ///
    /// For dim == 2 input, ipc::AABB always stores a 3-wide array whose z
    /// component is zero-initialized and never touched by the inflation (only
    /// the first `dim` components of the constructor argument are assigned) --
    /// so the z bound is an exact, uninflated 0.0. Replicate that exactly.
    ///
    /// @param vertices_t0 Positions at t0, column-major (dim * n).
    /// @param vertices_t1 Positions at t1, column-major (dim * n).
    /// @param n The number of vertices.
    /// @param dim The simulation dimension (2 or 3).
    /// @param inflation_radius The inflation radius.
    /// @param[out] box_min The box min corners (always 3 * n, row-major).
    /// @param[out] box_max The box max corners (always 3 * n, row-major).
    __global__ void build_vertex_boxes_kernel(
        const double* vertices_t0, // not __restrict__: may alias vertices_t1
        const double* vertices_t1, // (the static build passes one array twice)
        const int n,
        const int dim,
        const double inflation_radius,
        double* __restrict__ box_min,
        double* __restrict__ box_max)
    {
        const int i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= n) {
            return;
        }
#pragma unroll
        for (int k = 0; k < 3; ++k) {
            if (k < dim) {
                const double a = vertices_t0[k * n + i];
                const double b = vertices_t1[k * n + i];
                box_min[3 * i + k] = AABB::conservative_lower_bound(
                    fmin(a, b), inflation_radius);
                box_max[3 * i + k] = AABB::conservative_upper_bound(
                    fmax(a, b), inflation_radius);
            } else {
                box_min[3 * i + k] = 0.0;
                box_max[3 * i + k] = 0.0;
            }
        }
    }

    /// @brief One edge box as the union of its two vertex boxes, i.e. the
    /// AABB(aabb1, aabb2) constructor ipc::build_edge_boxes uses.
    __global__ void build_edge_boxes_kernel(
        const double* __restrict__ vbox_min,
        const double* __restrict__ vbox_max,
        const int32_t* __restrict__ edges, // 2 * n, row-major
        const int n,
        double* __restrict__ box_min,
        double* __restrict__ box_max)
    {
        const int i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= n) {
            return;
        }
        const int32_t e0 = edges[2 * i + 0];
        const int32_t e1 = edges[2 * i + 1];
#pragma unroll
        for (int k = 0; k < 3; ++k) {
            box_min[3 * i + k] =
                fmin(vbox_min[3 * e0 + k], vbox_min[3 * e1 + k]);
            box_max[3 * i + k] =
                fmax(vbox_max[3 * e0 + k], vbox_max[3 * e1 + k]);
        }
    }

    /// @brief One face box as the union of its three vertex boxes, i.e. the
    /// AABB(aabb1, aabb2, aabb3) constructor ipc::build_face_boxes uses.
    __global__ void build_face_boxes_kernel(
        const double* __restrict__ vbox_min,
        const double* __restrict__ vbox_max,
        const int32_t* __restrict__ faces, // 3 * n, row-major
        const int n,
        double* __restrict__ box_min,
        double* __restrict__ box_max)
    {
        const int i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= n) {
            return;
        }
        const int32_t f0 = faces[3 * i + 0];
        const int32_t f1 = faces[3 * i + 1];
        const int32_t f2 = faces[3 * i + 2];
#pragma unroll
        for (int k = 0; k < 3; ++k) {
            box_min[3 * i + k] = fmin(
                vbox_min[3 * f0 + k],
                fmin(vbox_min[3 * f1 + k], vbox_min[3 * f2 + k]));
            box_max[3 * i + k] = fmax(
                vbox_max[3 * f0 + k],
                fmax(vbox_max[3 * f1 + k], vbox_max[3 * f2 + k]));
        }
    }

    // -- Tree building ------------------------------------------------------

    /// @brief Compute one Morton code per box from its (normalized) center.
    /// Mirrors the compute_morton_codes block of ipc::LBVH::init_bvh.
    /// @param box_min The box min corners (3 * n, row-major).
    /// @param box_max The box max corners (3 * n, row-major).
    /// @param n The number of boxes.
    /// @param domain The Morton-normalization domain (device resident).
    /// @param dim The simulation dimension (2 or 3).
    /// @param[out] codes The Morton codes.
    /// @param[out] box_ids The box ids (the identity, to be sorted with codes).
    __global__ void compute_morton_codes_kernel(
        const double* __restrict__ box_min,
        const double* __restrict__ box_max,
        const int n,
        const Domain* __restrict__ domain,
        const int dim,
        uint64_t* __restrict__ codes,
        int32_t* __restrict__ box_ids)
    {
        const int i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= n) {
            return;
        }

        const Eigen::Array3d mesh_min(
            domain->min[0], domain->min[1], domain->min[2]);
        const Eigen::Array3d mesh_max(
            domain->max[0], domain->max[1], domain->max[2]);
        // The same derivation the CPU (ipc::LBVH::init_bvh) evaluates once per
        // build; IEEE division is deterministic, so evaluating it per thread
        // gives the identical reciprocal and so identical codes.
        const Eigen::Array3d mesh_width_inv =
            ipc::morton_domain_width_inv(mesh_min, mesh_max);

        const Eigen::Array3d center(
            0.5 * (box_min[3 * i + 0] + box_max[3 * i + 0]),
            0.5 * (box_min[3 * i + 1] + box_max[3 * i + 1]),
            0.5 * (box_min[3 * i + 2] + box_max[3 * i + 2]));

        codes[i] = ipc::morton_code(center, mesh_min, mesh_width_inv, dim);
        box_ids[i] = i;
    }

    /// @brief Single-pass bottom-up hierarchy + AABB build (Apetrei 2014).
    /// One thread per leaf, driving ipc::details::build_hierarchy_from_leaf()
    /// -- the same walk the CPU build runs -- with atomicAdd() and the fences
    /// around it standing in for the host's std::atomic arrival.
    /// @param box_min The box min corners (3 * n, row-major).
    /// @param box_max The box max corners (3 * n, row-major).
    /// @param sorted_codes The Morton codes in sorted order.
    /// @param sorted_box_ids The box ids in Morton-sorted order.
    /// @param n_leaves The number of leaves.
    /// @param[out] nodes The BVH nodes.
    /// @param[out] rightmost The per-node rightmost-leaf indices.
    /// @param[in,out] infos The per-node construction scratch (zeroed).
    /// @param[out] root_idx The root's index.
    __global__ void build_hierarchy_kernel(
        const double* __restrict__ box_min,
        const double* __restrict__ box_max,
        const uint64_t* __restrict__ sorted_codes,
        const int32_t* __restrict__ sorted_box_ids,
        const int n_leaves,
        ipc::LBVH::Node* __restrict__ nodes,
        int32_t* __restrict__ rightmost,
        DeviceConstructionInfo* __restrict__ infos,
        int* __restrict__ root_idx)
    {
        const int i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= n_leaves) {
            return;
        }

        const int32_t bid = sorted_box_ids[i];
        ipc::details::init_leaf_node(
            i, n_leaves, bid,
            Eigen::Array3d(
                box_min[3 * bid + 0], box_min[3 * bid + 1],
                box_min[3 * bid + 2]),
            Eigen::Array3d(
                box_max[3 * bid + 0], box_max[3 * bid + 1],
                box_max[3 * bid + 2]),
            nodes, rightmost);

        const int root = ipc::details::build_hierarchy_from_leaf(
            i, n_leaves, [sorted_codes](int k) { return sorted_codes[k]; },
            nodes, rightmost, infos,
            [](int& count) {
                // Release: publish this thread's child pointer, range endpoint
                // and leaf/subtree AABB before announcing arrival, so whoever
                // continues sees a complete child.
                __threadfence();
                const int previous = atomicAdd(&count, 1);
                if (previous != 0) {
                    // Acquire: this thread continues and reads the sibling's
                    // node, range endpoint and rightmost leaf. Those are
                    // ordinary loads, so without this fence they may be served
                    // from a stale L1 on another SM.
                    __threadfence();
                }
                return previous;
            });

        if (root >= 0) {
            *root_idx = root; // only one thread reaches the root
        }
    }

    /// @brief Swap the node and rightmost-leaf entries at index 0 and the root
    /// (single thread). See ipc::details::swap_root_to_zero(). Reads the root
    /// from device memory so the build never round-trips through the host.
    /// @param nodes The BVH nodes.
    /// @param rightmost The per-node rightmost-leaf indices.
    /// @param root The root's index; a no-op if it is already 0 (or unset).
    __global__ void swap_root_kernel(
        ipc::LBVH::Node* __restrict__ nodes,
        int32_t* __restrict__ rightmost,
        const int* __restrict__ root)
    {
        if (blockIdx.x == 0 && threadIdx.x == 0 && *root > 0) {
            ipc::details::swap_root_to_zero(nodes, rightmost, *root);
        }
    }

    /// @brief After the root swap, rewrite left pointers that referenced the
    /// old node 0 to its new location. See ipc::details::patch_left_pointer().
    /// @param nodes The BVH nodes.
    /// @param num_nodes The number of nodes.
    /// @param root The new location of the old node 0; a no-op if the root was
    /// already 0 (or unset).
    __global__ void patch_left_kernel(
        ipc::LBVH::Node* __restrict__ nodes,
        const int num_nodes,
        const int* __restrict__ root)
    {
        const int i = blockIdx.x * blockDim.x + threadIdx.x;
        const int r = *root;
        if (i >= num_nodes || r <= 0) {
            return;
        }
        ipc::details::patch_left_pointer(nodes[i], r);
    }

    /// @brief Build one BVH on the device from device-resident box corners.
    /// Mirrors ipc::LBVH::init_bvh. Fully asynchronous: nothing here waits on
    /// the device; the caller synchronizes once after all three trees.
    /// @param impl The pimpl whose scratch buffers to use.
    /// @param d_box_min The box min corners (3 * n, row-major, device).
    /// @param d_box_max The box max corners (3 * n, row-major, device).
    /// @param n The number of boxes (leaves).
    /// @param dim The simulation dimension (2 or 3).
    /// @param[out] bvh The BVH to build.
    /// @param[out] d_root Where to write the root's index (device); -1 if the
    /// build never reaches a root.
    void build_tree(
        LBVH::Impl& impl,
        const double* d_box_min,
        const double* d_box_max,
        const int n,
        const int dim,
        LBVH::Impl::DeviceBVH& bvh,
        int* d_root)
    {
        if (n == 0) {
            bvh.clear();
            return;
        }

        const size_t num_nodes = size_t(2) * n - 1;
        bvh.nodes.resize(num_nodes);
        bvh.rightmost_leaves.resize(num_nodes);

        for (auto& codes : impl.morton_codes) {
            codes.resize(n);
        }
        for (auto& ids : impl.box_ids) {
            ids.resize(n);
        }
        // Only the visitation counts need zeroing; a memset is the cheapest
        // way to do it.
        impl.construction_infos.resize(num_nodes);
        impl.construction_infos.zero();
        IPC_TOOLKIT_CUDA_CHECK(cudaMemsetAsync(d_root, 0xFF, sizeof(int)));

        compute_morton_codes_kernel<<<kernel_grid_size(n), KERNEL_BLOCK_SIZE>>>(
            d_box_min, d_box_max, n, impl.domain.data(), dim,
            impl.morton_codes[0].data(), impl.box_ids[0].data());
        IPC_TOOLKIT_CUDA_CHECK(cudaGetLastError());

        // Radix sort the (code, id) pairs by code.
        //
        // A radix sort is out of place: each pass reads one array and writes
        // another, so CUB is given two buffers per sequence (a DoubleBuffer)
        // and ping-pongs between them, pass after pass. The unsorted input is
        // in buffer [0], where the codes kernel wrote it; after the sort,
        // Current() says which of the two holds the result -- it depends on
        // the number of passes, so it must be read back rather than assumed.
        cub::DoubleBuffer<uint64_t> keys(
            impl.morton_codes[0].data(), impl.morton_codes[1].data());
        cub::DoubleBuffer<int32_t> values(
            impl.box_ids[0].data(), impl.box_ids[1].data());

        // CUB's two-phase convention: with a null storage pointer the call
        // only reports the scratch bytes it needs; the second call sorts.
        // Keeping that scratch in a persistent buffer is what makes this
        // allocation-free per build (thrust::sort_by_key would cudaMalloc and
        // cudaFree it every call).
        size_t temp_bytes = 0;
        IPC_TOOLKIT_CUDA_CHECK(
            cub::DeviceRadixSort::SortPairs(
                nullptr, temp_bytes, keys, values, n));
        impl.sort_temp.resize(temp_bytes);
        IPC_TOOLKIT_CUDA_CHECK(
            cub::DeviceRadixSort::SortPairs(
                impl.sort_temp.data(), temp_bytes, keys, values, n));

        // Current() is the sorted half of each ping-pong pair.
        build_hierarchy_kernel<<<kernel_grid_size(n), KERNEL_BLOCK_SIZE>>>(
            d_box_min, d_box_max, keys.Current(), values.Current(), n,
            bvh.nodes.data(), bvh.rightmost_leaves.data(),
            impl.construction_infos.data(), d_root);
        IPC_TOOLKIT_CUDA_CHECK(cudaGetLastError());

        swap_root_kernel<<<1, 1>>>(
            bvh.nodes.data(), bvh.rightmost_leaves.data(), d_root);
        IPC_TOOLKIT_CUDA_CHECK(cudaGetLastError());

        patch_left_kernel<<<kernel_grid_size(num_nodes), KERNEL_BLOCK_SIZE>>>(
            bvh.nodes.data(), static_cast<int>(num_nodes), d_root);
        IPC_TOOLKIT_CUDA_CHECK(cudaGetLastError());
    }

    /// @brief Reduce the Morton-normalization domain (min of mins, max of
    /// maxs) over the device-resident vertex box corners into impl.domain.
    /// The same quantity BroadPhase::compute_mesh_aabb() computes on the host
    /// from the same boxes, with the same seeds; min/max are order-independent,
    /// so the two agree bit-for-bit. Stays on the device: the codes kernel
    /// reads it directly, with no host round-trip.
    void compute_domain(LBVH::Impl& impl, const int n_vertices)
    {
        Domain init;
        for (int k = 0; k < 3; ++k) {
            init.min[k] = std::numeric_limits<double>::max();
            init.max[k] = std::numeric_limits<double>::lowest();
        }

        const auto domains = thrust::make_transform_iterator(
            thrust::counting_iterator<int>(0),
            MakeDomain { impl.vbox_min.data(), impl.vbox_max.data() });

        impl.domain.resize(1);
        size_t temp_bytes = 0;
        IPC_TOOLKIT_CUDA_CHECK(
            cub::DeviceReduce::Reduce(
                nullptr, temp_bytes, domains, impl.domain.data(), n_vertices,
                DomainReduce {}, init));
        impl.reduce_temp.resize(temp_bytes);
        IPC_TOOLKIT_CUDA_CHECK(
            cub::DeviceReduce::Reduce(
                impl.reduce_temp.data(), temp_bytes, domains,
                impl.domain.data(), n_vertices, DomainReduce {}, init));
    }

    /// @brief Upload vertex positions column-major, straight from the matrix
    /// when it is contiguous (no host transpose).
    void upload_vertices(
        Eigen::ConstRef<Eigen::MatrixXd> vertices, DeviceBuffer<double>& d)
    {
        const size_t n = static_cast<size_t>(vertices.size());
        if (vertices.innerStride() == 1
            && vertices.outerStride() == vertices.rows()) {
            d.upload(vertices.data(), n);
        } else {
            const Eigen::MatrixXd contiguous = vertices;
            d.upload(contiguous.data(), n);
        }
    }

    /// @brief Flatten an integer connectivity matrix (rowwise) to row-major
    /// 32-bit ids on the host, and upload the same array to the device.
    template <int Cols>
    void upload_connectivity(
        Eigen::ConstRef<Eigen::MatrixXi> M,
        std::vector<int32_t>& h,
        DeviceBuffer<int32_t>& d)
    {
        const size_t n = M.rows();
        h.resize(Cols * n);
        for (size_t i = 0; i < n; ++i) {
            for (int k = 0; k < Cols; ++k) {
                h[Cols * i + k] = static_cast<int32_t>(M(i, k));
            }
        }
        d.upload(h.data(), h.size());
    }

    /// @brief Copy a device BVH to the host. Bytewise: ipc::LBVH::Node holds
    /// only floats and ints, so a memcpy is its copy.
    void to_host(
        const LBVH::Impl::DeviceBVH& bvh,
        ipc::LBVH::Nodes& nodes,
        ipc::LBVH::RightmostLeaves& rightmost_leaves)
    {
        nodes.resize(bvh.nodes.size());
        rightmost_leaves.resize(bvh.rightmost_leaves.size());
        bvh.nodes.download(nodes.data());
        bvh.rightmost_leaves.download(rightmost_leaves.data());
    }

    /// @brief Given device-resident vertex boxes (impl.vbox_min/max), build the
    /// edge/face boxes and all three BVHs. Shared by every build() overload.
    /// @param impl The pimpl to fill (output).
    /// @param dim The simulation dimension (2 or 3).
    /// @param n_vertices The number of vertices.
    /// @param edges The mesh edges.
    /// @param faces The mesh faces.
    void build_from_vertex_boxes(
        LBVH::Impl& impl,
        const int dim,
        const int n_vertices,
        Eigen::ConstRef<Eigen::MatrixXi> edges,
        Eigen::ConstRef<Eigen::MatrixXi> faces)
    {
        assert(edges.size() == 0 || edges.cols() == 2);
        assert(faces.size() == 0 || faces.cols() == 3);

        const int n_edges = static_cast<int>(edges.rows());
        const int n_faces = static_cast<int>(faces.rows());

        upload_connectivity<2>(edges, impl.h_edges, impl.edges);
        upload_connectivity<3>(faces, impl.h_faces, impl.faces);

        // Build edge/face boxes on the device from the vertex boxes.
        impl.ebox_min.resize(3 * size_t(n_edges));
        impl.ebox_max.resize(3 * size_t(n_edges));
        if (n_edges > 0) {
            build_edge_boxes_kernel<<<
                kernel_grid_size(n_edges), KERNEL_BLOCK_SIZE>>>(
                impl.vbox_min.data(), impl.vbox_max.data(), impl.edges.data(),
                n_edges, impl.ebox_min.data(), impl.ebox_max.data());
            IPC_TOOLKIT_CUDA_CHECK(cudaGetLastError());
        }

        impl.fbox_min.resize(3 * size_t(n_faces));
        impl.fbox_max.resize(3 * size_t(n_faces));
        if (n_faces > 0) {
            build_face_boxes_kernel<<<
                kernel_grid_size(n_faces), KERNEL_BLOCK_SIZE>>>(
                impl.vbox_min.data(), impl.vbox_max.data(), impl.faces.data(),
                n_faces, impl.fbox_min.data(), impl.fbox_max.data());
            IPC_TOOLKIT_CUDA_CHECK(cudaGetLastError());
        }

        // The CPU normalizes all three BVHs by the vertex box domain.
        compute_domain(impl, n_vertices);

        impl.roots.resize(3);
        build_tree(
            impl, impl.vbox_min.data(), impl.vbox_max.data(), n_vertices, dim,
            impl.vertex_bvh, impl.roots.data() + 0);
        build_tree(
            impl, impl.ebox_min.data(), impl.ebox_max.data(), n_edges, dim,
            impl.edge_bvh, impl.roots.data() + 1);
        build_tree(
            impl, impl.fbox_min.data(), impl.fbox_max.data(), n_faces, dim,
            impl.face_bvh, impl.roots.data() + 2);

        // The one synchronization of the build: surfaces any kernel fault and
        // lets the roots be read back -- all three at once.
        IPC_TOOLKIT_CUDA_CHECK(cudaDeviceSynchronize());

        // A hierarchy build that never reaches its root leaves -1 behind. The
        // traversal always starts at node 0, which would then be an arbitrary
        // interior node, silently dropping every candidate outside its
        // subtree. Refuse to hand out such a tree.
        int roots[3];
        impl.roots.download(roots);
        const int n_leaves[3] = { n_vertices, n_edges, n_faces };
        for (int k = 0; k < 3; ++k) {
            if (n_leaves[k] > 0 && roots[k] < 0) {
                log_and_throw_error(
                    "ipc::cuda::LBVH: the device hierarchy build did not "
                    "reach a root (tree {} of 3, {} leaves); the BVH is "
                    "malformed",
                    k, n_leaves[k]);
            }
        }
    }

    // -- Traversal ----------------------------------------------------------

    /// @brief Load a primitive's vertex ids into registers: a vertex's id is
    /// itself; an edge's or a face's come from its connectivity row. Fully
    /// unrolled, so the array is only ever indexed by constants and stays in
    /// registers (a runtime-indexed local array would go to local memory).
    /// @tparam N The number of vertex ids per primitive (1, 2, or 3).
    /// @param prim The primitive id.
    /// @param conn The connectivity (N ids per primitive, row-major); unused
    /// and may be null for N == 1.
    /// @param[out] ids The primitive's vertex ids.
    template <int N>
    __device__ inline void load_vertex_ids(
        const int32_t prim, const int32_t* __restrict__ conn, int32_t (&ids)[N])
    {
        if constexpr (N == 1) {
            ids[0] = prim;
        } else {
            assert(conn != nullptr);
#pragma unroll
            for (int k = 0; k < N; ++k) {
                ids[k] = conn[N * prim + k];
            }
        }
    }

    /// @brief Append a (source_prim, target_prim) pair (post-swap) via an
    /// atomic counter. Writes only if the slot is within capacity; the counter
    /// still advances on overflow so the caller learns the required size. The
    /// counter and capacity are 64-bit so neither can wrap or truncate.
    template <bool swap_order>
    __device__ inline void emit_pair(
        const int32_t query_prim,
        const int32_t node_prim,
        int32_t* __restrict__ out_a,
        int32_t* __restrict__ out_b,
        unsigned long long* __restrict__ counter,
        const unsigned long long capacity)
    {
        int32_t a = query_prim, b = node_prim;
        if constexpr (swap_order) {
            const int32_t t = a;
            a = b;
            b = t;
        }
        const unsigned long long slot = atomicAdd(counter, 1ULL);
        if (slot < capacity) {
            out_a[slot] = a;
            out_b[slot] = b;
        }
    }

    /// @brief One thread per source leaf: descend the target BVH and append
    /// every AABB-overlapping, connectivity-passing (source_prim, target_prim)
    /// pair to the output arrays. The descent is ipc::details::traverse_lbvh()
    /// and the shared-vertex exclusion ipc::details::share_vertex(), both
    /// shared with the CPU ipc::LBVH. The remaining user vertex filter (if
    /// any) is applied on the host, so the final set matches the CPU.
    /// @tparam triangular Self-collision: skip subtrees left of the query.
    /// @tparam swap_order Emit (target_prim, source_prim) instead.
    /// @tparam SourceCount The vertex ids per source primitive (1, 2, or 3).
    /// @tparam TargetCount The vertex ids per target primitive (1, 2, or 3).
    /// @param source The BVH whose leaves are the queries.
    /// @param n_source_leaves The number of source leaves.
    /// @param source_leaf_offset The index of the source BVH's first leaf.
    /// @param target The BVH to descend.
    /// @param target_size The number of nodes in the target BVH.
    /// @param target_rightmost The target's per-node rightmost-leaf indices.
    /// @param source_conn The source connectivity (null for vertices).
    /// @param target_conn The target connectivity (null for vertices).
    /// @param[out] out_a The first ids of the emitted pairs.
    /// @param[out] out_b The second ids of the emitted pairs.
    /// @param[in,out] counter The emitted-pair counter.
    /// @param capacity The output arrays' capacity.
    template <
        bool triangular,
        bool swap_order,
        int SourceCount,
        int TargetCount>
    __global__ void traverse_kernel(
        const ipc::LBVH::Node* __restrict__ source,
        const int n_source_leaves,
        const int source_leaf_offset,
        const ipc::LBVH::Node* __restrict__ target,
        const int target_size,
        const int32_t* __restrict__ target_rightmost,
        const int32_t* __restrict__ source_conn,
        const int32_t* __restrict__ target_conn,
        int32_t* __restrict__ out_a,
        int32_t* __restrict__ out_b,
        unsigned long long* __restrict__ counter,
        const unsigned long long capacity)
    {
        const int s = blockIdx.x * blockDim.x + threadIdx.x;
        if (s >= n_source_leaves) {
            return;
        }
        const ipc::LBVH::Node query = source[source_leaf_offset + s];

        int32_t query_ids[SourceCount];
        load_vertex_ids<SourceCount>(
            query.primitive_id, source_conn, query_ids);

        ipc::details::traverse_lbvh<triangular>(
            s, target, target_size, target_rightmost,
            [&](const ipc::LBVH::Node& node) { return node.intersects(query); },
            [&](const ipc::LBVH::Node& leaf, const int /*leaf_idx*/,
                const bool /*intersects*/) {
                int32_t leaf_ids[TargetCount];
                load_vertex_ids<TargetCount>(
                    leaf.primitive_id, target_conn, leaf_ids);
                if (!ipc::details::share_vertex(query_ids, leaf_ids)) {
                    emit_pair<swap_order>(
                        query.primitive_id, leaf.primitive_id, out_a, out_b,
                        counter, capacity);
                }
            });
    }

    /// @brief How one candidate type is traversed: which BVH is the source
    /// (its leaves are the queries), which is the target (descended), how many
    /// vertex ids each primitive has, whether the pair is triangular (a BVH
    /// against itself) and whether it is emitted swapped. Each tuple is stated
    /// exactly once here and drives both the host-materializing and the
    /// device-view detect paths, so the two cannot disagree.
    template <typename Candidate> struct Traversal;

    template <> struct Traversal<VertexVertexCandidate> {
        static constexpr bool triangular = true;
        static constexpr bool swap_order = false;
        static constexpr int source_count = 1;
        static constexpr int target_count = 1;
        static const LBVH::Impl::DeviceBVH& source(const LBVH::Impl& impl)
        {
            return impl.vertex_bvh;
        }
        static const LBVH::Impl::DeviceBVH& target(const LBVH::Impl& impl)
        {
            return impl.vertex_bvh;
        }
        static const int32_t* source_conn(const LBVH::Impl&) { return nullptr; }
        static const int32_t* target_conn(const LBVH::Impl&) { return nullptr; }
        static LBVH::Impl::DeviceCandidates& buffer(LBVH::Impl& impl)
        {
            return impl.vv_candidates;
        }
    };

    // In 2D and for codimensional edge-vertex collisions there are more
    // vertices than edges, so iterate over the edges. Mirrors ipc::LBVH.
    template <> struct Traversal<EdgeVertexCandidate> {
        static constexpr bool triangular = false;
        static constexpr bool swap_order = false;
        static constexpr int source_count = 2;
        static constexpr int target_count = 1;
        static const LBVH::Impl::DeviceBVH& source(const LBVH::Impl& impl)
        {
            return impl.edge_bvh;
        }
        static const LBVH::Impl::DeviceBVH& target(const LBVH::Impl& impl)
        {
            return impl.vertex_bvh;
        }
        static const int32_t* source_conn(const LBVH::Impl& impl)
        {
            return impl.edges.data();
        }
        static const int32_t* target_conn(const LBVH::Impl&) { return nullptr; }
        static LBVH::Impl::DeviceCandidates& buffer(LBVH::Impl& impl)
        {
            return impl.ev_candidates;
        }
    };

    template <> struct Traversal<EdgeEdgeCandidate> {
        static constexpr bool triangular = true;
        static constexpr bool swap_order = false;
        static constexpr int source_count = 2;
        static constexpr int target_count = 2;
        static const LBVH::Impl::DeviceBVH& source(const LBVH::Impl& impl)
        {
            return impl.edge_bvh;
        }
        static const LBVH::Impl::DeviceBVH& target(const LBVH::Impl& impl)
        {
            return impl.edge_bvh;
        }
        static const int32_t* source_conn(const LBVH::Impl& impl)
        {
            return impl.edges.data();
        }
        static const int32_t* target_conn(const LBVH::Impl& impl)
        {
            return impl.edges.data();
        }
        static LBVH::Impl::DeviceCandidates& buffer(LBVH::Impl& impl)
        {
            return impl.ee_candidates;
        }
    };

    // The ratio vertices:faces is 1:2, so iterate over the vertices and query
    // the face BVH, swapping so the emitted pair is (face, vertex). Mirrors
    // ipc::LBVH.
    template <> struct Traversal<FaceVertexCandidate> {
        static constexpr bool triangular = false;
        static constexpr bool swap_order = true;
        static constexpr int source_count = 1;
        static constexpr int target_count = 3;
        static const LBVH::Impl::DeviceBVH& source(const LBVH::Impl& impl)
        {
            return impl.vertex_bvh;
        }
        static const LBVH::Impl::DeviceBVH& target(const LBVH::Impl& impl)
        {
            return impl.face_bvh;
        }
        static const int32_t* source_conn(const LBVH::Impl&) { return nullptr; }
        static const int32_t* target_conn(const LBVH::Impl& impl)
        {
            return impl.faces.data();
        }
        static LBVH::Impl::DeviceCandidates& buffer(LBVH::Impl& impl)
        {
            return impl.fv_candidates;
        }
    };

    // The ratio edges:faces is 3:2, so iterate over the faces and query the
    // edge BVH, swapping so the emitted pair is (edge, face). Mirrors
    // ipc::LBVH.
    template <> struct Traversal<EdgeFaceCandidate> {
        static constexpr bool triangular = false;
        static constexpr bool swap_order = true;
        static constexpr int source_count = 3;
        static constexpr int target_count = 2;
        static const LBVH::Impl::DeviceBVH& source(const LBVH::Impl& impl)
        {
            return impl.face_bvh;
        }
        static const LBVH::Impl::DeviceBVH& target(const LBVH::Impl& impl)
        {
            return impl.edge_bvh;
        }
        static const int32_t* source_conn(const LBVH::Impl& impl)
        {
            return impl.faces.data();
        }
        static const int32_t* target_conn(const LBVH::Impl& impl)
        {
            return impl.edges.data();
        }
        static LBVH::Impl::DeviceCandidates& buffer(LBVH::Impl& impl)
        {
            return impl.ef_candidates;
        }
    };

    template <> struct Traversal<FaceFaceCandidate> {
        static constexpr bool triangular = true;
        static constexpr bool swap_order = false;
        static constexpr int source_count = 3;
        static constexpr int target_count = 3;
        static const LBVH::Impl::DeviceBVH& source(const LBVH::Impl& impl)
        {
            return impl.face_bvh;
        }
        static const LBVH::Impl::DeviceBVH& target(const LBVH::Impl& impl)
        {
            return impl.face_bvh;
        }
        static const int32_t* source_conn(const LBVH::Impl& impl)
        {
            return impl.faces.data();
        }
        static const int32_t* target_conn(const LBVH::Impl& impl)
        {
            return impl.faces.data();
        }
        static LBVH::Impl::DeviceCandidates& buffer(LBVH::Impl& impl)
        {
            return impl.ff_candidates;
        }
    };

    /// @brief Run the device traversal for one candidate type, leaving the
    /// connectivity-filtered pairs device-resident in its buffer.
    ///
    /// The first pass is sized from the buffer's high-water mark (the largest
    /// count any earlier call on this object needed), or a guess of 8 pairs
    /// per source leaf, whichever is larger. The kernel counts every pair it
    /// finds but writes only those that fit, so if the first pass overflows
    /// the exact count is known and a second pass always fits: at most two
    /// passes, and only the first call -- or a call whose count exceeds every
    /// prior one -- pays for the second.
    ///
    /// @tparam Candidate The candidate type; see Traversal.
    /// @param impl The pimpl; the caller must hold impl.mutex.
    /// @return The number of candidate pairs emitted.
    template <typename Candidate> size_t run_traversal(LBVH::Impl& impl)
    {
        using T = Traversal<Candidate>;
        const LBVH::Impl::DeviceBVH& source = T::source(impl);
        const LBVH::Impl::DeviceBVH& target = T::target(impl);
        LBVH::Impl::DeviceCandidates& buf = T::buffer(impl);

        buf.clear();

        const int n_source_leaves = source.n_leaves();
        const int target_size = static_cast<int>(target.nodes.size());
        // A triangular traversal is a BVH against itself, and a lone primitive
        // cannot collide with itself.
        if (n_source_leaves == 0 || target_size == 0
            || (T::triangular && n_source_leaves < 2)) {
            return 0;
        }
        const int source_leaf_offset = n_source_leaves - 1;

        size_t capacity = std::max(
            { buf.a.capacity(), size_t(1024),
              size_t(8) * static_cast<size_t>(n_source_leaves) });
        impl.counter.resize(1);

        unsigned long long count = 0;
        for (int pass = 0; pass < 2; ++pass) {
            buf.a.reserve(capacity);
            buf.b.reserve(capacity);
            impl.counter.zero();

            traverse_kernel<
                T::triangular, T::swap_order, T::source_count, T::target_count>
                <<<kernel_grid_size(n_source_leaves), KERNEL_BLOCK_SIZE>>>(
                    source.nodes.data(), n_source_leaves, source_leaf_offset,
                    target.nodes.data(), target_size,
                    target.rightmost_leaves.data(), T::source_conn(impl),
                    T::target_conn(impl), buf.a.data(), buf.b.data(),
                    impl.counter.data(),
                    static_cast<unsigned long long>(capacity));
            IPC_TOOLKIT_CUDA_CHECK(cudaGetLastError());

            impl.counter.download(&count); // synchronizes
            if (count <= capacity) {
                break; // everything fit
            }

            logger().warn(
                "[ipc::cuda::LBVH] candidate count {} exceeded preallocated "
                "capacity {}; re-running with the exact size (this cost is "
                "amortized: later calls reuse the learned capacity)",
                count, capacity);
            capacity = static_cast<size_t>(count); // exact size now known
        }
        if (count > capacity) {
            // The count is a deterministic function of the two trees, so the
            // second pass must fit; anything else is a device fault.
            log_and_throw_error(
                "ipc::cuda::LBVH: the candidate count changed between two "
                "identical traversals ({} > capacity {})",
                count, capacity);
        }

        buf.count = static_cast<size_t>(count);
        buf.a.resize(buf.count);
        buf.b.resize(buf.count);
        return buf.count;
    }

    /// @brief Copy the device-resident candidate pairs to host Candidate
    /// objects. For the accept-all filter every pair is kept (the device set
    /// is already exact); otherwise the user vertex filter trims the
    /// connectivity-filtered superset.
    /// @param buf The device pairs.
    /// @param filter The user vertex filter.
    /// @param can_collide The full predicate applied when the filter is not accept-all.
    /// @param[out] out The materialized candidates (cleared first).
    template <typename Candidate, typename CanCollide>
    void materialize(
        const LBVH::Impl::DeviceCandidates& buf,
        const CollisionFilter& filter,
        const CanCollide& can_collide,
        std::vector<Candidate>& out)
    {
        out.clear();
        const size_t count = buf.count;
        if (count == 0) {
            return;
        }
        std::vector<int32_t> h_a(count), h_b(count);
        buf.a.download(h_a.data());
        buf.b.download(h_b.data());

        out.reserve(count);
        if (filter.accepts_all()) {
            for (size_t k = 0; k < count; ++k) {
                out.emplace_back(h_a[k], h_b[k]);
            }
        } else {
            for (size_t k = 0; k < count; ++k) {
                if (can_collide(h_a[k], h_b[k])) {
                    out.emplace_back(h_a[k], h_b[k]);
                }
            }
        }
    }

    /// @brief Traverse on the device, then materialize on the host.
    template <typename Candidate, typename CanCollide>
    void detect_host(
        LBVH::Impl& impl,
        const CollisionFilter& filter,
        const CanCollide& can_collide,
        std::vector<Candidate>& out)
    {
        // The detect_*() methods are const on the BroadPhase interface, and
        // ipc::LBVH's really are read-only, so a caller may legitimately run
        // two of them concurrently on one shared object. Here every one of
        // them writes shared device state -- the pair counter, and the
        // candidate buffer of its type -- so without this lock two detections
        // would race on the counter and on buffer reallocation. Held for the
        // materialize too, so another call cannot overwrite the buffer while
        // it is being copied to the host.
        const std::lock_guard<std::mutex> lock(impl.mutex);
        run_traversal<Candidate>(impl);
        materialize(
            Traversal<Candidate>::buffer(impl), filter, can_collide, out);
    }

    /// @brief Traverse on the device and return a view of the result.
    template <typename Candidate>
    LBVH::DeviceCandidateView detect_device(LBVH::Impl& impl)
    {
        // Same reason as detect_host(): const on the interface, but writes the
        // shared counter and this type's candidate buffer. The lock ends with
        // the call, so the returned view is only as safe as the caller's own
        // ordering of later detect_*() calls (see DeviceCandidateView).
        const std::lock_guard<std::mutex> lock(impl.mutex);
        const size_t count = run_traversal<Candidate>(impl);
        const LBVH::Impl::DeviceCandidates& buf =
            Traversal<Candidate>::buffer(impl);
        return LBVH::DeviceCandidateView { count ? buf.a.data() : nullptr,
                                           count ? buf.b.data() : nullptr,
                                           count };
    }

} // namespace

// ---------------------------------------------------------------------------

LBVH::LBVH() : ipc::BroadPhase(), m_impl(std::make_unique<Impl>()) { }

LBVH::~LBVH() = default;

// ipc::BroadPhase declares a destructor and so has no move operations: moving
// it as a whole would invoke its copy, which copies a std::function and may
// throw. Its members are moved individually instead, which cannot. The
// moved-from object is left in the cleared state (no Impl, dim 0, default
// filter); impl() re-seeds it on its next use.

LBVH::LBVH(LBVH&& other) noexcept
    : ipc::BroadPhase()
    , m_impl(std::move(other.m_impl))
{
    can_vertices_collide = std::move(other.can_vertices_collide);
    other.can_vertices_collide = CollisionFilter();
    dim = other.dim;
    other.dim = 0;
}

LBVH& LBVH::operator=(LBVH&& other) noexcept
{
    if (this != &other) {
        m_impl = std::move(other.m_impl);
        can_vertices_collide = std::move(other.can_vertices_collide);
        other.can_vertices_collide = CollisionFilter();
        dim = other.dim;
        other.dim = 0;
    }
    return *this;
}

LBVH::Impl& LBVH::impl() const
{
    if (!m_impl) {
        m_impl = std::make_unique<Impl>();
    }
    return *m_impl;
}

void LBVH::build(
    Eigen::ConstRef<Eigen::MatrixXd> vertices,
    Eigen::ConstRef<Eigen::MatrixXi> edges,
    Eigen::ConstRef<Eigen::MatrixXi> faces,
    const double inflation_radius)
{
    // A static box is a temporal one whose endpoints coincide. The dynamic
    // build notices the two refs alias and uploads the vertices only once.
    build(vertices, vertices, edges, faces, inflation_radius);
}

void LBVH::build(
    Eigen::ConstRef<Eigen::MatrixXd> vertices_t0,
    Eigen::ConstRef<Eigen::MatrixXd> vertices_t1,
    Eigen::ConstRef<Eigen::MatrixXi> edges,
    Eigen::ConstRef<Eigen::MatrixXi> faces,
    const double inflation_radius)
{
    assert(vertices_t0.rows() == vertices_t1.rows());
    assert(vertices_t0.cols() == vertices_t1.cols());
    assert(vertices_t0.rows() <= std::numeric_limits<int32_t>::max());

    clear();

    assert(vertices_t0.cols() == 2 || vertices_t0.cols() == 3);
    dim = static_cast<uint8_t>(vertices_t0.cols());

    const int n_vertices = static_cast<int>(vertices_t0.rows());
    if (n_vertices == 0) {
        return;
    }

    Impl& device = impl();
    upload_vertices(vertices_t0, device.vertices_t0);

    // The static build passes the same matrix twice; upload it once and point
    // the kernel's t1 at the t0 copy rather than paying a second transfer.
    const bool same_vertices = vertices_t0.data() == vertices_t1.data()
        && vertices_t0.outerStride() == vertices_t1.outerStride();
    if (!same_vertices) {
        upload_vertices(vertices_t1, device.vertices_t1);
    }
    const double* d_vertices_t1 =
        same_vertices ? device.vertices_t0.data() : device.vertices_t1.data();

    // Build vertex boxes on the device (always 3-wide storage).
    device.vbox_min.resize(3 * size_t(n_vertices));
    device.vbox_max.resize(3 * size_t(n_vertices));
    build_vertex_boxes_kernel<<<
        kernel_grid_size(n_vertices), KERNEL_BLOCK_SIZE>>>(
        device.vertices_t0.data(), d_vertices_t1, n_vertices, dim,
        inflation_radius, device.vbox_min.data(), device.vbox_max.data());
    IPC_TOOLKIT_CUDA_CHECK(cudaGetLastError());

    build_from_vertex_boxes(device, dim, n_vertices, edges, faces);
}

void LBVH::build(
    Eigen::ConstRef<Eigen::MatrixXi> edges,
    Eigen::ConstRef<Eigen::MatrixXi> faces)
{
    // BroadPhase::build(const AABBs&, edges, faces, dim) has cleared us and
    // filled vertex_boxes and dim.
    assert(dim == 2 || dim == 3);
    assert(vertex_boxes.size() <= std::numeric_limits<int32_t>::max());

    const int n_vertices = static_cast<int>(vertex_boxes.size());
    if (n_vertices == 0) {
        return;
    }

    // Upload the precomputed vertex boxes.
    std::vector<double> h_min(3 * size_t(n_vertices));
    std::vector<double> h_max(3 * size_t(n_vertices));
    for (int i = 0; i < n_vertices; ++i) {
        for (int k = 0; k < 3; ++k) {
            h_min[3 * size_t(i) + k] = vertex_boxes[i].min[k];
            h_max[3 * size_t(i) + k] = vertex_boxes[i].max[k];
        }
    }
    Impl& device = impl();
    device.vbox_min.upload(h_min.data(), h_min.size());
    device.vbox_max.upload(h_max.data(), h_max.size());

    build_from_vertex_boxes(device, dim, n_vertices, edges, faces);

    // As in ipc::LBVH: the host boxes are redundant once the trees exist.
    vertex_boxes.clear();
}

void LBVH::clear()
{
    ipc::BroadPhase::clear();
    if (m_impl) { // a moved-from object has nothing to clear
        m_impl->clear();
    }
}

// ---------------------------------------------------------------------------
// BroadPhase interface. Device BVH descent + device connectivity filter; the
// user vertex filter is applied on the host only when it is not accept-all.

void LBVH::detect_vertex_vertex_candidates(
    std::vector<VertexVertexCandidate>& candidates) const
{
    detect_host<VertexVertexCandidate>(
        impl(), can_vertices_collide,
        [this](size_t a, size_t b) { return can_vertices_collide(a, b); },
        candidates);
}

void LBVH::detect_edge_vertex_candidates(
    std::vector<EdgeVertexCandidate>& candidates) const
{
    detect_host<EdgeVertexCandidate>(
        impl(), can_vertices_collide,
        [this](size_t a, size_t b) { return can_edge_vertex_collide(a, b); },
        candidates);
}

void LBVH::detect_edge_edge_candidates(
    std::vector<EdgeEdgeCandidate>& candidates) const
{
    detect_host<EdgeEdgeCandidate>(
        impl(), can_vertices_collide,
        [this](size_t a, size_t b) { return can_edges_collide(a, b); },
        candidates);
}

void LBVH::detect_face_vertex_candidates(
    std::vector<FaceVertexCandidate>& candidates) const
{
    detect_host<FaceVertexCandidate>(
        impl(), can_vertices_collide,
        [this](size_t a, size_t b) { return can_face_vertex_collide(a, b); },
        candidates);
}

void LBVH::detect_edge_face_candidates(
    std::vector<EdgeFaceCandidate>& candidates) const
{
    detect_host<EdgeFaceCandidate>(
        impl(), can_vertices_collide,
        [this](size_t a, size_t b) { return can_edge_face_collide(a, b); },
        candidates);
}

void LBVH::detect_face_face_candidates(
    std::vector<FaceFaceCandidate>& candidates) const
{
    detect_host<FaceFaceCandidate>(
        impl(), can_vertices_collide,
        [this](size_t a, size_t b) { return can_faces_collide(a, b); },
        candidates);
}

// ---------------------------------------------------------------------------
// Device-resident candidate accessors.

LBVH::DeviceCandidateView LBVH::detect_vertex_vertex_candidates_device() const
{
    return detect_device<VertexVertexCandidate>(impl());
}

LBVH::DeviceCandidateView LBVH::detect_edge_vertex_candidates_device() const
{
    return detect_device<EdgeVertexCandidate>(impl());
}

LBVH::DeviceCandidateView LBVH::detect_edge_edge_candidates_device() const
{
    return detect_device<EdgeEdgeCandidate>(impl());
}

LBVH::DeviceCandidateView LBVH::detect_face_vertex_candidates_device() const
{
    return detect_device<FaceVertexCandidate>(impl());
}

LBVH::DeviceCandidateView LBVH::detect_edge_face_candidates_device() const
{
    return detect_device<EdgeFaceCandidate>(impl());
}

LBVH::DeviceCandidateView LBVH::detect_face_face_candidates_device() const
{
    return detect_device<FaceFaceCandidate>(impl());
}

// ---------------------------------------------------------------------------
// Host-side can_*_collide filters (mesh connectivity + user vertex filter).
// Mirror ipc::LBVH's overrides, backed by the host connectivity copies. The
// ids arrive from device memory, so the bounds asserts are the tripwire that
// localizes a host/device desync to the traversal.

bool LBVH::can_edge_vertex_collide(size_t ei, size_t vi) const
{
    const std::vector<int32_t>& edges = impl().h_edges;
    assert(2 * ei + 1 < edges.size());

    return ipc::details::can_edge_vertex_collide(
        edges[2 * ei], edges[2 * ei + 1], vi, can_vertices_collide);
}

bool LBVH::can_edges_collide(size_t eai, size_t ebi) const
{
    const std::vector<int32_t>& edges = impl().h_edges;
    assert(2 * eai + 1 < edges.size());
    assert(2 * ebi + 1 < edges.size());

    return ipc::details::can_edges_collide(
        edges[2 * eai], edges[2 * eai + 1], edges[2 * ebi], edges[2 * ebi + 1],
        can_vertices_collide);
}

bool LBVH::can_face_vertex_collide(size_t fi, size_t vi) const
{
    const std::vector<int32_t>& faces = impl().h_faces;
    assert(3 * fi + 2 < faces.size());

    return ipc::details::can_face_vertex_collide(
        faces[3 * fi], faces[3 * fi + 1], faces[3 * fi + 2], vi,
        can_vertices_collide);
}

bool LBVH::can_edge_face_collide(size_t ei, size_t fi) const
{
    const std::vector<int32_t>& edges = impl().h_edges;
    const std::vector<int32_t>& faces = impl().h_faces;
    assert(2 * ei + 1 < edges.size());
    assert(3 * fi + 2 < faces.size());

    return ipc::details::can_edge_face_collide(
        edges[2 * ei], edges[2 * ei + 1], faces[3 * fi], faces[3 * fi + 1],
        faces[3 * fi + 2], can_vertices_collide);
}

bool LBVH::can_faces_collide(size_t fai, size_t fbi) const
{
    const std::vector<int32_t>& faces = impl().h_faces;
    assert(3 * fai + 2 < faces.size());
    assert(3 * fbi + 2 < faces.size());

    return ipc::details::can_faces_collide(
        faces[3 * fai], faces[3 * fai + 1], faces[3 * fai + 2], faces[3 * fbi],
        faces[3 * fbi + 1], faces[3 * fbi + 2], can_vertices_collide);
}

size_t LBVH::num_vertex_nodes() const { return impl().vertex_bvh.nodes.size(); }

size_t LBVH::num_edge_nodes() const { return impl().edge_bvh.nodes.size(); }

size_t LBVH::num_face_nodes() const { return impl().face_bvh.nodes.size(); }

// ---------------------------------------------------------------------------
// Debug / validation.

void LBVH::vertex_nodes_to_host(
    ipc::LBVH::Nodes& nodes, ipc::LBVH::RightmostLeaves& rightmost_leaves) const
{
    to_host(impl().vertex_bvh, nodes, rightmost_leaves);
}

void LBVH::edge_nodes_to_host(
    ipc::LBVH::Nodes& nodes, ipc::LBVH::RightmostLeaves& rightmost_leaves) const
{
    to_host(impl().edge_bvh, nodes, rightmost_leaves);
}

void LBVH::face_nodes_to_host(
    ipc::LBVH::Nodes& nodes, ipc::LBVH::RightmostLeaves& rightmost_leaves) const
{
    to_host(impl().face_bvh, nodes, rightmost_leaves);
}

} // namespace ipc::cuda

#endif // IPC_TOOLKIT_WITH_CUDA
