#include "lbvh.hpp"

#ifdef IPC_TOOLKIT_WITH_CUDA

#include <ipc/broad_phase/cuda/lbvh_impl.cuh>
#include <ipc/math/morton.hpp>
#include <ipc/utils/cuda/device_utils.cuh>
#include <ipc/utils/logger.hpp>

#include <thrust/copy.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/sort.h>
#include <thrust/transform_reduce.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <vector>

namespace ipc::cuda {

namespace {

    // Eigen::Array3d is passed to kernels by value, so it must be exactly three
    // packed doubles (no vectorization padding) to have a stable layout.
    static_assert(
        sizeof(Eigen::Array3d) == 24,
        "Eigen::Array3d must be 24 bytes (3 packed doubles)");

    /// @brief Per-internal-node scratch used by the bottom-up build. The device
    /// analog of ipc::LBVH::ConstructionInfo, kept separate on purpose: that
    /// struct's visitation_count is a std::atomic<int>, which cannot be used
    /// here (atomicAdd needs an int*, and std::atomic is non-copyable so it
    /// cannot be a thrust::device_vector element). A plain int suffices because
    /// atomicAdd provides the atomicity the CPU gets from std::atomic.
    struct DeviceConstructionInfo {
        int left_range;
        int right_range;
        int visitation_count;
    };

    /// @brief Min/max domain accumulator for the Morton-normalization reduction.
    struct Domain {
        double mn[3];
        double mx[3];
    };

    struct DomainReduce {
        __host__ __device__ Domain
        operator()(const Domain& a, const Domain& b) const
        {
            Domain r;
#pragma unroll
            for (int k = 0; k < 3; ++k) {
                r.mn[k] = fmin(a.mn[k], b.mn[k]);
                r.mx[k] = fmax(a.mx[k], b.mx[k]);
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
                d.mn[k] = box_min[3 * i + k];
                d.mx[k] = box_max[3 * i + k];
            }
            return d;
        }
    };

    // -- Box building -------------------------------------------------------
    // Matches ipc::build_*_boxes + AABB::conservative_inflation exactly: the
    // double bounds are nudged outward with nextafter so the box is
    // conservative. (The leaf nodes later apply a second float-nextafter in
    // build_hierarchy_kernel, matching assign_inflated_aabb.)
    //
    // For dim == 2 input, ipc::AABB always stores a 3-wide array whose z
    // component is zero-initialized and never touched by conservative_inflation
    // (only the first `dim` components of the constructor argument are
    // assigned) -- so the z bound is an exact, uninflated 0.0, not
    // nextafter(0 +/- inflation_radius, ...). Replicate that exactly: for
    // k >= dim, write a hard 0.0 instead of inflating.

    __global__ void build_vertex_boxes_static_kernel(
        const double* __restrict__ vertices, // dim * n, row-major
        const int n,
        const int dim,
        const double inflation_radius,
        double* __restrict__ box_min, // always 3 * n, row-major
        double* __restrict__ box_max)
    {
        const int i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= n) {
            return;
        }
#pragma unroll
        for (int k = 0; k < 3; ++k) {
            if (k < dim) {
                const double v = vertices[dim * i + k];
                box_min[3 * i + k] = nextafter(v - inflation_radius, -INFINITY);
                box_max[3 * i + k] = nextafter(v + inflation_radius, INFINITY);
            } else {
                box_min[3 * i + k] = 0.0;
                box_max[3 * i + k] = 0.0;
            }
        }
    }

    __global__ void build_vertex_boxes_dynamic_kernel(
        const double* __restrict__ vertices_t0, // dim * n, row-major
        const double* __restrict__ vertices_t1, // dim * n, row-major
        const int n,
        const int dim,
        const double inflation_radius,
        double* __restrict__ box_min, // always 3 * n, row-major
        double* __restrict__ box_max)
    {
        const int i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= n) {
            return;
        }
#pragma unroll
        for (int k = 0; k < 3; ++k) {
            if (k < dim) {
                const double a = vertices_t0[dim * i + k];
                const double b = vertices_t1[dim * i + k];
                // union of the two inflated point boxes; nextafter is
                // monotonic so min(nextafter(a),nextafter(b)) ==
                // nextafter(min(a,b)).
                box_min[3 * i + k] =
                    nextafter(fmin(a, b) - inflation_radius, -INFINITY);
                box_max[3 * i + k] =
                    nextafter(fmax(a, b) + inflation_radius, INFINITY);
            } else {
                box_min[3 * i + k] = 0.0;
                box_max[3 * i + k] = 0.0;
            }
        }
    }

    __global__ void build_edge_boxes_kernel(
        const double* __restrict__ vbox_min,
        const double* __restrict__ vbox_max,
        const index_t* __restrict__ edges, // 2 * n, row-major
        const int n,
        double* __restrict__ box_min,
        double* __restrict__ box_max)
    {
        const int i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= n) {
            return;
        }
        const index_t e0 = edges[2 * i + 0];
        const index_t e1 = edges[2 * i + 1];
#pragma unroll
        for (int k = 0; k < 3; ++k) {
            box_min[3 * i + k] =
                fmin(vbox_min[3 * e0 + k], vbox_min[3 * e1 + k]);
            box_max[3 * i + k] =
                fmax(vbox_max[3 * e0 + k], vbox_max[3 * e1 + k]);
        }
    }

    __global__ void build_face_boxes_kernel(
        const double* __restrict__ vbox_min,
        const double* __restrict__ vbox_max,
        const index_t* __restrict__ faces, // 3 * n, row-major
        const int n,
        double* __restrict__ box_min,
        double* __restrict__ box_max)
    {
        const int i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= n) {
            return;
        }
        const index_t f0 = faces[3 * i + 0];
        const index_t f1 = faces[3 * i + 1];
        const index_t f2 = faces[3 * i + 2];
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

    /// @brief Number of common leading bits between Morton codes at sorted
    /// positions i and j (device port of the CPU delta()). Duplicate codes fall
    /// back to the CLZ of the index XOR (offset by 32 so it sorts after any
    /// code-level difference).
    /// @param sorted_codes The Morton codes in sorted order.
    /// @param n The number of codes.
    /// @param i The first sorted position.
    /// @param code_i The code at position i (passed to avoid a redundant look-up).
    /// @param j The second sorted position.
    /// @return The common-prefix length, or -1 when j is out of bounds.
    __device__ inline int delta_device(
        const uint64_t* __restrict__ sorted_codes,
        const int n,
        const int i,
        const uint64_t code_i,
        const int j)
    {
        if (j < 0 || j >= n) {
            return -1;
        }
        const uint64_t code_j = sorted_codes[j];
        if (code_i == code_j) {
            return 32 + __clz(i ^ j);
        }
        return __clzll(static_cast<long long>(code_i ^ code_j));
    }

    /// @brief Compute one Morton code per box from its (normalized) center.
    /// Mirrors the compute_morton_codes block of ipc::LBVH::init_bvh.
    __global__ void compute_morton_codes_kernel(
        const double* __restrict__ box_min, // 3 * n, row-major
        const double* __restrict__ box_max, // 3 * n, row-major
        const int n,
        const Eigen::Array3d mesh_min,
        const Eigen::Array3d mesh_width_inv,
        const int dim,
        uint64_t* __restrict__ codes,
        index_t* __restrict__ box_ids)
    {
        const int i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= n) {
            return;
        }

        const double cx = 0.5 * (box_min[3 * i + 0] + box_max[3 * i + 0]);
        const double cy = 0.5 * (box_min[3 * i + 1] + box_max[3 * i + 1]);
        const double cz = 0.5 * (box_min[3 * i + 2] + box_max[3 * i + 2]);

        // (center - mesh_min) * mesh_width_inv -- the reciprocal is
        // precomputed once per build (see compute_domain) and multiplied here
        // instead of dividing per box, matching the CPU (ipc::LBVH::init_bvh)
        // bit-for-bit.
        const double mx = (cx - mesh_min.x()) * mesh_width_inv.x();
        const double my = (cy - mesh_min.y()) * mesh_width_inv.y();
        const double mz = (cz - mesh_min.z()) * mesh_width_inv.z();

        codes[i] = (dim == 2) ? morton_2D(mx, my) : morton_3D(mx, my, mz);
        box_ids[i] = i;
    }

    /// @brief Single-pass bottom-up hierarchy + AABB build (Apetrei 2014).
    /// One thread per leaf. Direct port of the build_hierarchy_and_boxes block
    /// of ipc::LBVH::init_bvh, with atomicAdd + __threadfence replacing the
    /// std::atomic arrival gate.
    __global__ void build_hierarchy_kernel(
        const double* __restrict__ box_min,
        const double* __restrict__ box_max,
        const uint64_t* __restrict__ sorted_codes,
        const index_t* __restrict__ sorted_box_ids,
        const int N_LEAVES,
        ipc::LBVH::Node* __restrict__ nodes,
        int32_t* __restrict__ rightmost,
        DeviceConstructionInfo* __restrict__ infos,
        int* __restrict__ root_idx)
    {
        const int i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= N_LEAVES) {
            return;
        }

        const int LEAF_OFFSET = N_LEAVES - 1;

        // --- Initialize leaf node ---
        {
            const index_t bid = sorted_box_ids[i];
            ipc::LBVH::Node leaf;
#pragma unroll
            for (int k = 0; k < 3; ++k) {
                // Round the float AABB out (matches assign_inflated_aabb).
                leaf.aabb_min[k] = nextafterf(
                    static_cast<float>(box_min[3 * bid + k]), -INFINITY);
                leaf.aabb_max[k] = nextafterf(
                    static_cast<float>(box_max[3 * bid + k]), INFINITY);
            }
            leaf.primitive_id = static_cast<int32_t>(bid);
            leaf.is_inner_marker = 0;
            nodes[LEAF_OFFSET + i] = leaf;
            // A leaf's rightmost leaf is itself.
            rightmost[LEAF_OFFSET + i] = i;
        }

        // Single-node tree: the leaf is the root; no internal nodes to build.
        if (N_LEAVES == 1) {
            if (i == 0) {
                *root_idx = 0;
            }
            return;
        }

        // --- Bottom-up walk (Apetrei 2014, Fig. 2) ---
        int left_key = i;
        int right_key = i;
        int current_node = LEAF_OFFSET + i;

        while (true) {
            // Choose parent (see the CPU comment in ipc::LBVH::init_bvh).
            const bool is_child_a = (left_key == 0)
                || (right_key != N_LEAVES - 1
                    && delta_device(
                           sorted_codes, N_LEAVES, right_key,
                           sorted_codes[right_key], right_key + 1)
                        > delta_device(
                            sorted_codes, N_LEAVES, left_key - 1,
                            sorted_codes[left_key - 1], left_key));
            const int parent = is_child_a ? right_key : left_key - 1;

            // Write the child pointer + range onto the parent.
            if (is_child_a) {
                nodes[parent].left = current_node;
                infos[parent].left_range = left_key;
            } else {
                nodes[parent].right = current_node;
                infos[parent].right_range = right_key;
            }

            // Publish this child's node data and range to all threads before
            // signaling arrival, so the second thread reads consistent state.
            __threadfence();

            // Atomic arrival gate: first thread stops; second proceeds knowing
            // both children are complete.
            if (atomicAdd(&infos[parent].visitation_count, 1) == 0) {
                break; // first thread to arrive -> finished
            }

            // Second thread: compute the parent AABB union and rightmost leaf.
            const ipc::LBVH::Node& child_a = nodes[nodes[parent].left];
            const ipc::LBVH::Node& child_b = nodes[nodes[parent].right];
            nodes[parent].aabb_min = child_a.aabb_min.min(child_b.aabb_min);
            nodes[parent].aabb_max = child_a.aabb_max.max(child_b.aabb_max);
            rightmost[parent] = ::max(
                rightmost[nodes[parent].left], rightmost[nodes[parent].right]);

            // Reconstruct the parent's full key range and continue upward.
            left_key = infos[parent].left_range;
            right_key = infos[parent].right_range;
            current_node = parent;

            if (left_key == 0 && right_key == N_LEAVES - 1) {
                // Only one thread reaches the root.
                *root_idx = current_node;
                break;
            }
        }
    }

    /// @brief Swap the node and rightmost-leaf entries at indices 0 and root
    /// (runs on a single thread).
    /// @param nodes The BVH nodes.
    /// @param rightmost The per-node rightmost-leaf indices.
    /// @param root The index to swap with index 0.
    __global__ void swap_root_kernel(
        ipc::LBVH::Node* __restrict__ nodes,
        int32_t* __restrict__ rightmost,
        const int root)
    {
        if (blockIdx.x == 0 && threadIdx.x == 0) {
            const ipc::LBVH::Node tmp = nodes[0];
            nodes[0] = nodes[root];
            nodes[root] = tmp;
            const int32_t t = rightmost[0];
            rightmost[0] = rightmost[root];
            rightmost[root] = t;
        }
    }

    /// @brief After the root swap, rewrite left pointers that referenced the
    /// old node 0 to its new location. See the CPU swap_root_to_zero comment:
    /// the old node 0 was only ever a left child, so only .left needs patching.
    /// is_inner_marker aliases .right and is nonzero iff internal.
    /// @param nodes The BVH nodes.
    /// @param num_nodes The number of nodes.
    /// @param root The new location of the old node 0.
    __global__ void patch_left_kernel(
        ipc::LBVH::Node* __restrict__ nodes,
        const int num_nodes,
        const int root)
    {
        const int i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= num_nodes) {
            return;
        }
        if (nodes[i].is_inner_marker != 0 && nodes[i].left == 0) {
            nodes[i].left = root;
        }
    }

    /// @brief Build one BVH on the device from device-resident box corners.
    /// Mirrors ipc::LBVH::init_bvh; the output BVH is resized and filled in
    /// place.
    /// @param d_box_min The box min corners (3 * n, row-major, device).
    /// @param d_box_max The box max corners (3 * n, row-major, device).
    /// @param n The number of boxes (leaves).
    /// @param mesh_min The Morton-normalization domain minimum.
    /// @param mesh_width_inv The reciprocal of the Morton-normalization domain
    /// extent (precomputed once per build; see compute_domain).
    /// @param dim The simulation dimension (2 or 3).
    /// @param bvh The BVH to build (output).
    void build_tree(
        const double* d_box_min,
        const double* d_box_max,
        const int n,
        const Eigen::Array3d& mesh_min,
        const Eigen::Array3d& mesh_width_inv,
        const int dim,
        LBVH::Impl::DeviceBVH& bvh)
    {
        bvh.n_leaves = n;
        if (n == 0) {
            bvh.nodes.clear();
            bvh.rightmost_leaves.clear();
            return;
        }

        const size_t num_nodes = size_t(2) * n - 1;
        bvh.nodes.resize(num_nodes);
        bvh.rightmost_leaves.resize(num_nodes);

        thrust::device_vector<uint64_t> morton_codes(n);
        thrust::device_vector<index_t> box_ids(n);
        // Value-initialized to zero => visitation_count starts at 0.
        thrust::device_vector<DeviceConstructionInfo> infos(num_nodes);
        thrust::device_vector<int> d_root(1, -1);

        compute_morton_codes_kernel<<<kernel_grid_size(n), KERNEL_BLOCK_SIZE>>>(
            d_box_min, d_box_max, n, mesh_min, mesh_width_inv, dim,
            thrust::raw_pointer_cast(morton_codes.data()),
            thrust::raw_pointer_cast(box_ids.data()));
        IPC_TOOLKIT_CUDA_CHECK(cudaGetLastError());

        thrust::sort_by_key(
            morton_codes.begin(), morton_codes.end(), box_ids.begin());

        build_hierarchy_kernel<<<kernel_grid_size(n), KERNEL_BLOCK_SIZE>>>(
            d_box_min, d_box_max, thrust::raw_pointer_cast(morton_codes.data()),
            thrust::raw_pointer_cast(box_ids.data()), n,
            thrust::raw_pointer_cast(bvh.nodes.data()),
            thrust::raw_pointer_cast(bvh.rightmost_leaves.data()),
            thrust::raw_pointer_cast(infos.data()),
            thrust::raw_pointer_cast(d_root.data()));
        IPC_TOOLKIT_CUDA_CHECK(cudaGetLastError());

        const int root = d_root[0]; // device->host read
        if (root > 0) {
            swap_root_kernel<<<1, 1>>>(
                thrust::raw_pointer_cast(bvh.nodes.data()),
                thrust::raw_pointer_cast(bvh.rightmost_leaves.data()), root);
            IPC_TOOLKIT_CUDA_CHECK(cudaGetLastError());

            patch_left_kernel<<<
                kernel_grid_size(num_nodes), KERNEL_BLOCK_SIZE>>>(
                thrust::raw_pointer_cast(bvh.nodes.data()),
                static_cast<int>(num_nodes), root);
            IPC_TOOLKIT_CUDA_CHECK(cudaGetLastError());
        }
    }

    /// @brief Compute the Morton-normalization domain (min of mins, max of
    /// maxs) over device-resident vertex box corners, and its reciprocal
    /// extent. The reciprocal is computed once here (per build) and multiplied
    /// per box in compute_morton_codes_kernel instead of dividing per box,
    /// matching the CPU (ipc::LBVH::init_bvh) bit-for-bit.
    void compute_domain(
        const thrust::device_vector<double>& vbox_min,
        const thrust::device_vector<double>& vbox_max,
        const int n_vertices,
        Eigen::Array3d& mesh_min,
        Eigen::Array3d& mesh_width_inv)
    {
        Domain init;
        for (int k = 0; k < 3; ++k) {
            init.mn[k] = std::numeric_limits<double>::max();
            init.mx[k] = std::numeric_limits<double>::lowest();
        }
        const Domain dom = thrust::transform_reduce(
            thrust::counting_iterator<int>(0),
            thrust::counting_iterator<int>(n_vertices),
            MakeDomain { thrust::raw_pointer_cast(vbox_min.data()),
                         thrust::raw_pointer_cast(vbox_max.data()) },
            init, DomainReduce {});

        mesh_min = Eigen::Array3d(dom.mn[0], dom.mn[1], dom.mn[2]);
        const Eigen::Array3d mesh_width(
            dom.mx[0] - dom.mn[0], dom.mx[1] - dom.mn[1],
            dom.mx[2] - dom.mn[2]);
        mesh_width_inv = 1.0 / mesh_width;
    }

    // Upload an integer connectivity matrix (rowwise) as a flat row-major
    // index_t device array.
    template <int Cols>
    thrust::device_vector<index_t>
    upload_connectivity(Eigen::ConstRef<Eigen::MatrixXi> M)
    {
        const size_t n = M.rows();
        std::vector<index_t> h(Cols * n);
        for (size_t i = 0; i < n; ++i) {
            for (int k = 0; k < Cols; ++k) {
                h[Cols * i + k] = static_cast<index_t>(M(i, k));
            }
        }
        return thrust::device_vector<index_t>(h);
    }

    void to_host(
        const LBVH::Impl::DeviceBVH& bvh,
        ipc::LBVH::Nodes& nodes,
        ipc::LBVH::RightmostLeaves& rightmost_leaves)
    {
        nodes.resize(bvh.nodes.size());
        rightmost_leaves.resize(bvh.rightmost_leaves.size());
        thrust::copy(bvh.nodes.begin(), bvh.nodes.end(), nodes.begin());
        thrust::copy(
            bvh.rightmost_leaves.begin(), bvh.rightmost_leaves.end(),
            rightmost_leaves.begin());
    }

    /// @brief Given device-resident vertex boxes, build the edge/face boxes and
    /// all three BVHs. Shared by every build() overload.
    /// @param impl The pimpl to fill (output).
    /// @param dim The simulation dimension (2 or 3).
    /// @param vbox_min The vertex box min corners (3 * n_vertices, device).
    /// @param vbox_max The vertex box max corners (3 * n_vertices, device).
    /// @param n_vertices The number of vertices.
    /// @param edges The mesh edges.
    /// @param faces The mesh faces.
    void build_from_vertex_boxes(
        LBVH::Impl& impl,
        const int dim,
        const thrust::device_vector<double>& vbox_min,
        const thrust::device_vector<double>& vbox_max,
        const int n_vertices,
        Eigen::ConstRef<Eigen::MatrixXi> edges,
        Eigen::ConstRef<Eigen::MatrixXi> faces)
    {
        assert(edges.size() == 0 || edges.cols() == 2);
        assert(faces.size() == 0 || faces.cols() == 3);

        const int n_edges = static_cast<int>(edges.rows());
        const int n_faces = static_cast<int>(faces.rows());

        // Upload connectivity to the device, and keep a host copy for the
        // host-side can_*_collide filters.
        impl.edges = upload_connectivity<2>(edges);
        impl.faces = upload_connectivity<3>(faces);

        impl.h_edge_vertex_ids.resize(n_edges);
        for (int i = 0; i < n_edges; ++i) {
            impl.h_edge_vertex_ids[i] = { { static_cast<index_t>(edges(i, 0)),
                                            static_cast<index_t>(
                                                edges(i, 1)) } };
        }
        impl.h_face_vertex_ids.resize(n_faces);
        for (int i = 0; i < n_faces; ++i) {
            impl.h_face_vertex_ids[i] = { { static_cast<index_t>(faces(i, 0)),
                                            static_cast<index_t>(faces(i, 1)),
                                            static_cast<index_t>(
                                                faces(i, 2)) } };
        }

        // Build edge/face boxes on the device from the vertex boxes.
        thrust::device_vector<double> ebox_min(3 * size_t(n_edges));
        thrust::device_vector<double> ebox_max(3 * size_t(n_edges));
        if (n_edges > 0) {
            build_edge_boxes_kernel<<<
                kernel_grid_size(n_edges), KERNEL_BLOCK_SIZE>>>(
                thrust::raw_pointer_cast(vbox_min.data()),
                thrust::raw_pointer_cast(vbox_max.data()),
                thrust::raw_pointer_cast(impl.edges.data()), n_edges,
                thrust::raw_pointer_cast(ebox_min.data()),
                thrust::raw_pointer_cast(ebox_max.data()));
            IPC_TOOLKIT_CUDA_CHECK(cudaGetLastError());
        }

        thrust::device_vector<double> fbox_min(3 * size_t(n_faces));
        thrust::device_vector<double> fbox_max(3 * size_t(n_faces));
        if (n_faces > 0) {
            build_face_boxes_kernel<<<
                kernel_grid_size(n_faces), KERNEL_BLOCK_SIZE>>>(
                thrust::raw_pointer_cast(vbox_min.data()),
                thrust::raw_pointer_cast(vbox_max.data()),
                thrust::raw_pointer_cast(impl.faces.data()), n_faces,
                thrust::raw_pointer_cast(fbox_min.data()),
                thrust::raw_pointer_cast(fbox_max.data()));
            IPC_TOOLKIT_CUDA_CHECK(cudaGetLastError());
        }

        // The CPU normalizes all three BVHs by the vertex box domain.
        Eigen::Array3d mesh_min, mesh_width_inv;
        compute_domain(
            vbox_min, vbox_max, n_vertices, mesh_min, mesh_width_inv);

        build_tree(
            thrust::raw_pointer_cast(vbox_min.data()),
            thrust::raw_pointer_cast(vbox_max.data()), n_vertices, mesh_min,
            mesh_width_inv, dim, impl.vertex_bvh);
        build_tree(
            thrust::raw_pointer_cast(ebox_min.data()),
            thrust::raw_pointer_cast(ebox_max.data()), n_edges, mesh_min,
            mesh_width_inv, dim, impl.edge_bvh);
        build_tree(
            thrust::raw_pointer_cast(fbox_min.data()),
            thrust::raw_pointer_cast(fbox_max.data()), n_faces, mesh_min,
            mesh_width_inv, dim, impl.face_bvh);

        IPC_TOOLKIT_CUDA_CHECK(cudaDeviceSynchronize());
    }

    // -- Traversal ----------------------------------------------------------

    __device__ inline bool
    aabb_intersects(const ipc::LBVH::Node& a, const ipc::LBVH::Node& b)
    {
        return a.aabb_min[0] <= b.aabb_max[0] && b.aabb_min[0] <= a.aabb_max[0]
            && a.aabb_min[1] <= b.aabb_max[1] && b.aabb_min[1] <= a.aabb_max[1]
            && a.aabb_min[2] <= b.aabb_max[2] && b.aabb_min[2] <= a.aabb_max[2];
    }

    /// @brief Whether two primitives share a vertex id (the device connectivity
    /// filter). A vertex primitive's id set is {itself}; an edge's is its 2
    /// endpoints; a face's is its 3 vertices. This is exactly the
    /// shared-endpoint exclusion in ipc::LBVH's can_*_collide (for
    /// vertex-vertex it reduces to p_a == p_b).
    /// @param p_a The first primitive id.
    /// @param conn_a The first primitive's connectivity, or null for a vertex.
    /// @param count_a The number of vertex ids per first primitive (1, 2, or 3).
    /// @param p_b The second primitive id.
    /// @param conn_b The second primitive's connectivity, or null for a vertex.
    /// @param count_b The number of vertex ids per second primitive.
    /// @return Whether the two primitives share any vertex id.
    __device__ inline bool prim_shares_vertex(
        const int p_a,
        const index_t* __restrict__ conn_a,
        const int count_a,
        const int p_b,
        const index_t* __restrict__ conn_b,
        const int count_b)
    {
        index_t ids_a[3];
        index_t ids_b[3];
        if (conn_a == nullptr) {
            ids_a[0] = p_a;
        } else {
            for (int k = 0; k < count_a; ++k) {
                ids_a[k] = conn_a[count_a * p_a + k];
            }
        }
        if (conn_b == nullptr) {
            ids_b[0] = p_b;
        } else {
            for (int k = 0; k < count_b; ++k) {
                ids_b[k] = conn_b[count_b * p_b + k];
            }
        }
        for (int i = 0; i < count_a; ++i) {
            for (int j = 0; j < count_b; ++j) {
                if (ids_a[i] == ids_b[j]) {
                    return true;
                }
            }
        }
        return false;
    }

    /// @brief Append a (source_prim, target_prim) pair (post-swap) via an
    /// atomic counter. Writes only if the slot is within capacity; the counter
    /// still advances on overflow so the caller learns the required size.
    template <bool swap_order>
    __device__ inline void emit_pair(
        const int query_prim,
        const int node_prim,
        int32_t* __restrict__ out_a,
        int32_t* __restrict__ out_b,
        int* __restrict__ counter,
        const int capacity)
    {
        int a = query_prim, b = node_prim;
        if constexpr (swap_order) {
            const int t = a;
            a = b;
            b = t;
        }
        const int slot = atomicAdd(counter, 1);
        if (slot < capacity) {
            out_a[slot] = a;
            out_b[slot] = b;
        }
    }

    /// @brief One thread per source leaf: descend the target BVH and append
    /// every AABB-overlapping, connectivity-passing (source_prim, target_prim)
    /// pair to the output arrays. Descent is a direct port of traverse_lbvh()
    /// in lbvh.cpp (scalar path); the connectivity (shared-vertex) exclusion is
    /// applied here on the device. The remaining user vertex filter (if any) is
    /// applied on the host, so the final set matches the CPU ipc::LBVH.
    /// @tparam triangular Self-collision: skip subtrees fully left of the query.
    /// @tparam swap_order Emit (target_prim, source_prim) instead.
    template <bool triangular, bool swap_order>
    __global__ void traverse_kernel(
        const ipc::LBVH::Node* __restrict__ source,
        const int n_source_leaves,
        const int source_leaf_offset,
        const ipc::LBVH::Node* __restrict__ target,
        const int target_size,
        const int32_t* __restrict__ target_rightmost,
        const index_t* __restrict__ source_conn, // null for vertex primitives
        const int source_count,                  // ids per source primitive
        const index_t* __restrict__ target_conn, // null for vertex primitives
        const int target_count,                  // ids per target primitive
        int32_t* __restrict__ out_a,
        int32_t* __restrict__ out_b,
        int* __restrict__ counter,
        const int capacity)
    {
        const int s = blockIdx.x * blockDim.x + threadIdx.x;
        if (s >= n_source_leaves) {
            return;
        }
        const ipc::LBVH::Node query = source[source_leaf_offset + s];
        const int query_leaf_idx = s;

        constexpr int MAX_STACK_SIZE = 64;
        int stack[MAX_STACK_SIZE];
        int stack_ptr = 0;
        stack[stack_ptr++] = ipc::LBVH::Node::INVALID_POINTER; // 0

        int node_idx = 0; // root
        do {
            const ipc::LBVH::Node& node = target[node_idx];

            if (target_size == 1) { // single node (only root, which is a leaf)
                if constexpr (triangular) {
                    break; // no self-collision with a single primitive
                }
                if (aabb_intersects(node, query)
                    && !prim_shares_vertex(
                        query.primitive_id, source_conn, source_count,
                        node.primitive_id, target_conn, target_count)) {
                    emit_pair<swap_order>(
                        query.primitive_id, node.primitive_id, out_a, out_b,
                        counter, capacity);
                }
                break;
            }

            const ipc::LBVH::Node& child_l = target[node.left];
            const ipc::LBVH::Node& child_r = target[node.right];
            bool intersects_l = aabb_intersects(child_l, query);
            bool intersects_r = aabb_intersects(child_r, query);

            // Skip subtrees fully on the query's left (triangular only).
            if constexpr (triangular) {
                if (intersects_l
                    && target_rightmost[node.left] <= query_leaf_idx) {
                    intersects_l = false;
                }
                if (intersects_r
                    && target_rightmost[node.right] <= query_leaf_idx) {
                    intersects_r = false;
                }
            }

            // is_inner_marker aliases .right; it is 0 iff the node is a leaf.
            const bool l_leaf = (child_l.is_inner_marker == 0);
            const bool r_leaf = (child_r.is_inner_marker == 0);

            if (intersects_l && l_leaf
                && !prim_shares_vertex(
                    query.primitive_id, source_conn, source_count,
                    child_l.primitive_id, target_conn, target_count)) {
                emit_pair<swap_order>(
                    query.primitive_id, child_l.primitive_id, out_a, out_b,
                    counter, capacity);
            }
            if (intersects_r && r_leaf
                && !prim_shares_vertex(
                    query.primitive_id, source_conn, source_count,
                    child_r.primitive_id, target_conn, target_count)) {
                emit_pair<swap_order>(
                    query.primitive_id, child_r.primitive_id, out_a, out_b,
                    counter, capacity);
            }

            const bool traverse_l = intersects_l && !l_leaf;
            const bool traverse_r = intersects_r && !r_leaf;

            if (!traverse_l && !traverse_r) {
                node_idx = stack[--stack_ptr];
            } else {
                node_idx = traverse_l ? node.left : node.right;
                if (traverse_l && traverse_r) {
                    stack[stack_ptr++] = node.right;
                }
            }
        } while (node_idx != ipc::LBVH::Node::INVALID_POINTER);
    }

    /// @brief Run the device traversal of the target BVH by the source leaves,
    /// leaving the connectivity-filtered candidate pairs device-resident in
    /// buf.a/buf.b (resized to the exact count). The initial buffer size is
    /// seeded from buf.predicted_capacity (the largest count ever observed for
    /// this type on this object), so only the first call -- or a call whose
    /// count exceeds every prior call -- pays the overflow-and-retry cost;
    /// every other call fits on the first pass.
    /// @tparam triangular Self-collision: skip subtrees fully left of the query.
    /// @tparam swap_order Emit (target_prim, source_prim) instead.
    /// @param source The BVH whose leaves are the queries.
    /// @param target The BVH to descend.
    /// @param source_conn The source primitives' connectivity (null for vertices).
    /// @param source_count The vertex ids per source primitive (1, 2, or 3).
    /// @param target_conn The target primitives' connectivity (null for vertices).
    /// @param target_count The vertex ids per target primitive (1, 2, or 3).
    /// @param buf The output candidate buffer and capacity hint (in/out).
    /// @return The number of candidate pairs emitted.
    template <bool triangular, bool swap_order>
    size_t run_traversal(
        const LBVH::Impl::DeviceBVH& source,
        const LBVH::Impl::DeviceBVH& target,
        const index_t* source_conn,
        const int source_count,
        const index_t* target_conn,
        const int target_count,
        LBVH::Impl::DeviceCandidates& buf)
    {
        const int n_source_leaves = source.n_leaves;
        const int target_size = static_cast<int>(target.nodes.size());
        if (n_source_leaves == 0 || target_size == 0) {
            buf.a.clear();
            buf.b.clear();
            return 0;
        }
        const int source_leaf_offset = n_source_leaves - 1;

        size_t capacity = std::max(
            buf.predicted_capacity,
            static_cast<size_t>(std::max(1024, 8 * n_source_leaves)));
        thrust::device_vector<int> d_counter(1);

        int count = 0;
        while (true) {
            buf.a.resize(capacity);
            buf.b.resize(capacity);
            d_counter[0] = 0;

            traverse_kernel<triangular, swap_order>
                <<<kernel_grid_size(n_source_leaves), KERNEL_BLOCK_SIZE>>>(
                    thrust::raw_pointer_cast(source.nodes.data()),
                    n_source_leaves, source_leaf_offset,
                    thrust::raw_pointer_cast(target.nodes.data()), target_size,
                    thrust::raw_pointer_cast(target.rightmost_leaves.data()),
                    source_conn, source_count, target_conn, target_count,
                    thrust::raw_pointer_cast(buf.a.data()),
                    thrust::raw_pointer_cast(buf.b.data()),
                    thrust::raw_pointer_cast(d_counter.data()),
                    static_cast<int>(capacity));
            IPC_TOOLKIT_CUDA_CHECK(cudaGetLastError());

            count = d_counter[0]; // device->host read (also synchronizes)
            if (static_cast<size_t>(count) <= capacity) {
                break; // everything fit
            }

            logger().warn(
                "[ipc::cuda::LBVH] candidate count {} exceeded preallocated "
                "capacity {}; re-running with the exact size (this cost is "
                "amortized: later calls reuse the learned capacity)",
                count, capacity);
            capacity = static_cast<size_t>(count); // exact size now known
        }

        buf.predicted_capacity = std::max(buf.predicted_capacity, capacity);
        buf.a.resize(count); // shrink to the exact candidate count (keeps data)
        buf.b.resize(count);
        return static_cast<size_t>(count);
    }

    /// @brief Copy the device-resident candidate pairs to host Candidate
    /// objects. For the accept-all filter every pair is kept (the device set is
    /// already exact); otherwise the user vertex filter trims the
    /// connectivity-filtered superset.
    /// @param d_a The first ids of each candidate pair (device).
    /// @param d_b The second ids of each candidate pair (device).
    /// @param count The number of candidate pairs.
    /// @param accepts_all Whether the user vertex filter accepts every pair.
    /// @param can_collide The predicate applied when accepts_all is false.
    /// @param out The materialized candidates (appended to).
    template <typename Candidate>
    void materialize(
        const thrust::device_vector<int32_t>& d_a,
        const thrust::device_vector<int32_t>& d_b,
        const size_t count,
        const bool accepts_all,
        const std::function<bool(size_t, size_t)>& can_collide,
        std::vector<Candidate>& out)
    {
        if (count == 0) {
            return;
        }
        std::vector<int32_t> h_a(count), h_b(count);
        thrust::copy(d_a.begin(), d_a.begin() + count, h_a.begin());
        thrust::copy(d_b.begin(), d_b.begin() + count, h_b.begin());

        out.reserve(out.size() + count);
        if (accepts_all) {
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

} // namespace

LBVH::LBVH() : ipc::BroadPhase(), m_impl(std::make_unique<Impl>()) { }

LBVH::~LBVH() = default;

LBVH::LBVH(LBVH&&) noexcept = default;
LBVH& LBVH::operator=(LBVH&&) noexcept = default;

const LBVH::Impl& LBVH::impl() const { return *m_impl; }

void LBVH::build(
    Eigen::ConstRef<Eigen::MatrixXd> vertices,
    Eigen::ConstRef<Eigen::MatrixXi> edges,
    Eigen::ConstRef<Eigen::MatrixXi> faces,
    const double inflation_radius)
{
    clear();

    assert(vertices.cols() == 2 || vertices.cols() == 3);
    dim = static_cast<uint8_t>(vertices.cols());

    const int n_vertices = static_cast<int>(vertices.rows());
    if (n_vertices == 0) {
        return;
    }

    // Upload vertices as a flat row-major array (dim components per vertex;
    // no padding -- the box kernel below fills the unused z for 2D input).
    std::vector<double> h_verts(size_t(dim) * size_t(n_vertices));
    for (int i = 0; i < n_vertices; ++i) {
        for (int k = 0; k < dim; ++k) {
            h_verts[size_t(dim) * size_t(i) + k] = vertices(i, k);
        }
    }
    const thrust::device_vector<double> d_verts(h_verts);

    // Build vertex boxes on the device (always 3-wide storage).
    thrust::device_vector<double> vbox_min(3 * size_t(n_vertices));
    thrust::device_vector<double> vbox_max(3 * size_t(n_vertices));
    build_vertex_boxes_static_kernel<<<
        kernel_grid_size(n_vertices), KERNEL_BLOCK_SIZE>>>(
        thrust::raw_pointer_cast(d_verts.data()), n_vertices, dim,
        inflation_radius, thrust::raw_pointer_cast(vbox_min.data()),
        thrust::raw_pointer_cast(vbox_max.data()));
    IPC_TOOLKIT_CUDA_CHECK(cudaGetLastError());

    build_from_vertex_boxes(
        *m_impl, dim, vbox_min, vbox_max, n_vertices, edges, faces);
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

    clear();

    assert(vertices_t0.cols() == 2 || vertices_t0.cols() == 3);
    dim = static_cast<uint8_t>(vertices_t0.cols());

    const int n_vertices = static_cast<int>(vertices_t0.rows());
    if (n_vertices == 0) {
        return;
    }

    std::vector<double> h_v0(size_t(dim) * size_t(n_vertices));
    std::vector<double> h_v1(size_t(dim) * size_t(n_vertices));
    for (int i = 0; i < n_vertices; ++i) {
        for (int k = 0; k < dim; ++k) {
            h_v0[size_t(dim) * size_t(i) + k] = vertices_t0(i, k);
            h_v1[size_t(dim) * size_t(i) + k] = vertices_t1(i, k);
        }
    }
    const thrust::device_vector<double> d_v0(h_v0);
    const thrust::device_vector<double> d_v1(h_v1);

    thrust::device_vector<double> vbox_min(3 * size_t(n_vertices));
    thrust::device_vector<double> vbox_max(3 * size_t(n_vertices));
    build_vertex_boxes_dynamic_kernel<<<
        kernel_grid_size(n_vertices), KERNEL_BLOCK_SIZE>>>(
        thrust::raw_pointer_cast(d_v0.data()),
        thrust::raw_pointer_cast(d_v1.data()), n_vertices, dim,
        inflation_radius, thrust::raw_pointer_cast(vbox_min.data()),
        thrust::raw_pointer_cast(vbox_max.data()));
    IPC_TOOLKIT_CUDA_CHECK(cudaGetLastError());

    build_from_vertex_boxes(
        *m_impl, dim, vbox_min, vbox_max, n_vertices, edges, faces);
}

void LBVH::build(
    const AABBs& vertex_boxes,
    Eigen::ConstRef<Eigen::MatrixXi> edges,
    Eigen::ConstRef<Eigen::MatrixXi> faces,
    const uint8_t _dim)
{
    clear();

    assert(_dim == 2 || _dim == 3);
    dim = _dim;

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
    const thrust::device_vector<double> vbox_min(h_min);
    const thrust::device_vector<double> vbox_max(h_max);

    build_from_vertex_boxes(
        *m_impl, dim, vbox_min, vbox_max, n_vertices, edges, faces);
}

void LBVH::clear()
{
    ipc::BroadPhase::clear();
    if (m_impl) {
        m_impl->clear();
    }
}

// ---------------------------------------------------------------------------
// BroadPhase interface. Device BVH descent + device connectivity filter; the
// user vertex filter is applied on the host only when it is not accept-all.

namespace {
    // Raw device pointer to a connectivity array, or nullptr if empty (a
    // vertex primitive has no connectivity array).
    const index_t* conn_ptr(const thrust::device_vector<index_t>& v)
    {
        return v.empty() ? nullptr : thrust::raw_pointer_cast(v.data());
    }

    // Fill buf with the device connectivity-filtered candidate pairs, then
    // materialize them (host) into out, trimming with can_collide when the user
    // filter is not accept-all.
    template <typename Candidate, bool triangular, bool swap_order>
    void detect_host(
        const LBVH::Impl::DeviceBVH& source,
        const LBVH::Impl::DeviceBVH& target,
        const index_t* source_conn,
        const int source_count,
        const index_t* target_conn,
        const int target_count,
        LBVH::Impl::DeviceCandidates& buf,
        const bool accepts_all,
        const std::function<bool(size_t, size_t)>& can_collide,
        std::vector<Candidate>& out)
    {
        const size_t count = run_traversal<triangular, swap_order>(
            source, target, source_conn, source_count, target_conn,
            target_count, buf);
        materialize<Candidate>(
            buf.a, buf.b, count, accepts_all, can_collide, out);
    }

    // Fill buf on the device and return a view of it.
    template <bool triangular, bool swap_order>
    LBVH::DeviceCandidateView detect_device(
        const LBVH::Impl::DeviceBVH& source,
        const LBVH::Impl::DeviceBVH& target,
        const index_t* source_conn,
        const int source_count,
        const index_t* target_conn,
        const int target_count,
        LBVH::Impl::DeviceCandidates& buf)
    {
        const size_t count = run_traversal<triangular, swap_order>(
            source, target, source_conn, source_count, target_conn,
            target_count, buf);
        return LBVH::DeviceCandidateView {
            count ? thrust::raw_pointer_cast(buf.a.data()) : nullptr,
            count ? thrust::raw_pointer_cast(buf.b.data()) : nullptr, count
        };
    }
} // namespace

void LBVH::detect_vertex_vertex_candidates(
    std::vector<VertexVertexCandidate>& candidates) const
{
    if (m_impl->vertex_bvh.n_leaves <= 1) {
        return; // need at least 2 vertices for a collision
    }
    detect_host<VertexVertexCandidate, /*triangular=*/true, /*swap=*/false>(
        m_impl->vertex_bvh, m_impl->vertex_bvh, nullptr, 1, nullptr, 1,
        m_impl->vv_candidates, can_vertices_collide.accepts_all(),
        [this](size_t a, size_t b) { return can_vertices_collide(a, b); },
        candidates);
}

void LBVH::detect_edge_vertex_candidates(
    std::vector<EdgeVertexCandidate>& candidates) const
{
    if (m_impl->edge_bvh.n_leaves == 0 || m_impl->vertex_bvh.n_leaves == 0) {
        return;
    }
    detect_host<EdgeVertexCandidate, /*triangular=*/false, /*swap=*/false>(
        m_impl->edge_bvh, m_impl->vertex_bvh, conn_ptr(m_impl->edges), 2,
        nullptr, 1, m_impl->ev_candidates, can_vertices_collide.accepts_all(),
        [this](size_t a, size_t b) { return can_edge_vertex_collide(a, b); },
        candidates);
}

void LBVH::detect_edge_edge_candidates(
    std::vector<EdgeEdgeCandidate>& candidates) const
{
    if (m_impl->edge_bvh.n_leaves <= 1) {
        return; // need at least 2 edges for a collision
    }
    detect_host<EdgeEdgeCandidate, /*triangular=*/true, /*swap=*/false>(
        m_impl->edge_bvh, m_impl->edge_bvh, conn_ptr(m_impl->edges), 2,
        conn_ptr(m_impl->edges), 2, m_impl->ee_candidates,
        can_vertices_collide.accepts_all(),
        [this](size_t a, size_t b) { return can_edges_collide(a, b); },
        candidates);
}

void LBVH::detect_face_vertex_candidates(
    std::vector<FaceVertexCandidate>& candidates) const
{
    if (m_impl->face_bvh.n_leaves == 0 || m_impl->vertex_bvh.n_leaves == 0) {
        return;
    }
    // Iterate over the vertices (source) and query the face BVH (target),
    // swapping so the emitted pair is (face, vertex). Mirrors ipc::LBVH.
    detect_host<FaceVertexCandidate, /*triangular=*/false, /*swap=*/true>(
        m_impl->vertex_bvh, m_impl->face_bvh, nullptr, 1,
        conn_ptr(m_impl->faces), 3, m_impl->fv_candidates,
        can_vertices_collide.accepts_all(),
        [this](size_t a, size_t b) { return can_face_vertex_collide(a, b); },
        candidates);
}

void LBVH::detect_edge_face_candidates(
    std::vector<EdgeFaceCandidate>& candidates) const
{
    if (m_impl->edge_bvh.n_leaves == 0 || m_impl->face_bvh.n_leaves == 0) {
        return;
    }
    // Iterate over the faces (source) and query the edge BVH (target),
    // swapping so the emitted pair is (edge, face). Mirrors ipc::LBVH.
    detect_host<EdgeFaceCandidate, /*triangular=*/false, /*swap=*/true>(
        m_impl->face_bvh, m_impl->edge_bvh, conn_ptr(m_impl->faces), 3,
        conn_ptr(m_impl->edges), 2, m_impl->ef_candidates,
        can_vertices_collide.accepts_all(),
        [this](size_t a, size_t b) { return can_edge_face_collide(a, b); },
        candidates);
}

void LBVH::detect_face_face_candidates(
    std::vector<FaceFaceCandidate>& candidates) const
{
    if (m_impl->face_bvh.n_leaves <= 1) {
        return; // need at least 2 faces for a collision
    }
    detect_host<FaceFaceCandidate, /*triangular=*/true, /*swap=*/false>(
        m_impl->face_bvh, m_impl->face_bvh, conn_ptr(m_impl->faces), 3,
        conn_ptr(m_impl->faces), 3, m_impl->ff_candidates,
        can_vertices_collide.accepts_all(),
        [this](size_t a, size_t b) { return can_faces_collide(a, b); },
        candidates);
}

// ---------------------------------------------------------------------------
// Device-resident candidate accessors. Run the traversal and return a view of
// the connectivity-filtered pairs left on the device (valid until the next
// call on the same type or clear()). For the accept-all filter this is the
// exact candidate set; otherwise it is a superset the caller must trim with
// the user vertex filter.

LBVH::DeviceCandidateView LBVH::detect_vertex_vertex_candidates_device() const
{
    if (m_impl->vertex_bvh.n_leaves <= 1) {
        m_impl->vv_candidates.clear();
        return {};
    }
    return detect_device</*triangular=*/true, /*swap=*/false>(
        m_impl->vertex_bvh, m_impl->vertex_bvh, nullptr, 1, nullptr, 1,
        m_impl->vv_candidates);
}

LBVH::DeviceCandidateView LBVH::detect_edge_vertex_candidates_device() const
{
    if (m_impl->edge_bvh.n_leaves == 0 || m_impl->vertex_bvh.n_leaves == 0) {
        m_impl->ev_candidates.clear();
        return {};
    }
    return detect_device</*triangular=*/false, /*swap=*/false>(
        m_impl->edge_bvh, m_impl->vertex_bvh, conn_ptr(m_impl->edges), 2,
        nullptr, 1, m_impl->ev_candidates);
}

LBVH::DeviceCandidateView LBVH::detect_edge_edge_candidates_device() const
{
    if (m_impl->edge_bvh.n_leaves <= 1) {
        m_impl->ee_candidates.clear();
        return {};
    }
    return detect_device</*triangular=*/true, /*swap=*/false>(
        m_impl->edge_bvh, m_impl->edge_bvh, conn_ptr(m_impl->edges), 2,
        conn_ptr(m_impl->edges), 2, m_impl->ee_candidates);
}

LBVH::DeviceCandidateView LBVH::detect_face_vertex_candidates_device() const
{
    if (m_impl->face_bvh.n_leaves == 0 || m_impl->vertex_bvh.n_leaves == 0) {
        m_impl->fv_candidates.clear();
        return {};
    }
    return detect_device</*triangular=*/false, /*swap=*/true>(
        m_impl->vertex_bvh, m_impl->face_bvh, nullptr, 1,
        conn_ptr(m_impl->faces), 3, m_impl->fv_candidates);
}

LBVH::DeviceCandidateView LBVH::detect_edge_face_candidates_device() const
{
    if (m_impl->edge_bvh.n_leaves == 0 || m_impl->face_bvh.n_leaves == 0) {
        m_impl->ef_candidates.clear();
        return {};
    }
    return detect_device</*triangular=*/false, /*swap=*/true>(
        m_impl->face_bvh, m_impl->edge_bvh, conn_ptr(m_impl->faces), 3,
        conn_ptr(m_impl->edges), 2, m_impl->ef_candidates);
}

LBVH::DeviceCandidateView LBVH::detect_face_face_candidates_device() const
{
    if (m_impl->face_bvh.n_leaves <= 1) {
        m_impl->ff_candidates.clear();
        return {};
    }
    return detect_device</*triangular=*/true, /*swap=*/false>(
        m_impl->face_bvh, m_impl->face_bvh, conn_ptr(m_impl->faces), 3,
        conn_ptr(m_impl->faces), 3, m_impl->ff_candidates);
}

// ---------------------------------------------------------------------------
// Host-side can_*_collide filters (mesh connectivity + user vertex filter).
// Mirror ipc::LBVH's overrides, backed by the host connectivity copies.

bool LBVH::can_edge_vertex_collide(size_t ei, size_t vi) const
{
    const auto& [e0i, e1i] = m_impl->h_edge_vertex_ids[ei];
    return vi != e0i && vi != e1i
        && (can_vertices_collide(vi, e0i) || can_vertices_collide(vi, e1i));
}

bool LBVH::can_edges_collide(size_t eai, size_t ebi) const
{
    const auto& [ea0i, ea1i] = m_impl->h_edge_vertex_ids[eai];
    const auto& [eb0i, eb1i] = m_impl->h_edge_vertex_ids[ebi];

    const bool share_endpoint =
        ea0i == eb0i || ea0i == eb1i || ea1i == eb0i || ea1i == eb1i;

    return !share_endpoint
        && (can_vertices_collide(ea0i, eb0i) || can_vertices_collide(ea0i, eb1i)
            || can_vertices_collide(ea1i, eb0i)
            || can_vertices_collide(ea1i, eb1i));
}

bool LBVH::can_face_vertex_collide(size_t fi, size_t vi) const
{
    const auto& [f0i, f1i, f2i] = m_impl->h_face_vertex_ids[fi];
    return vi != f0i && vi != f1i && vi != f2i
        && (can_vertices_collide(vi, f0i) || can_vertices_collide(vi, f1i)
            || can_vertices_collide(vi, f2i));
}

bool LBVH::can_edge_face_collide(size_t ei, size_t fi) const
{
    const auto& [e0i, e1i] = m_impl->h_edge_vertex_ids[ei];
    const auto& [f0i, f1i, f2i] = m_impl->h_face_vertex_ids[fi];

    const bool share_endpoint = e0i == f0i || e0i == f1i || e0i == f2i
        || e1i == f0i || e1i == f1i || e1i == f2i;

    return !share_endpoint
        && (can_vertices_collide(e0i, f0i) || can_vertices_collide(e0i, f1i)
            || can_vertices_collide(e0i, f2i) || can_vertices_collide(e1i, f0i)
            || can_vertices_collide(e1i, f1i)
            || can_vertices_collide(e1i, f2i));
}

bool LBVH::can_faces_collide(size_t fai, size_t fbi) const
{
    const auto& [fa0i, fa1i, fa2i] = m_impl->h_face_vertex_ids[fai];
    const auto& [fb0i, fb1i, fb2i] = m_impl->h_face_vertex_ids[fbi];

    const bool share_endpoint = fa0i == fb0i || fa0i == fb1i || fa0i == fb2i
        || fa1i == fb0i || fa1i == fb1i || fa1i == fb2i || fa2i == fb0i
        || fa2i == fb1i || fa2i == fb2i;

    return !share_endpoint
        && (can_vertices_collide(fa0i, fb0i) || can_vertices_collide(fa0i, fb1i)
            || can_vertices_collide(fa0i, fb2i)
            || can_vertices_collide(fa1i, fb0i)
            || can_vertices_collide(fa1i, fb1i)
            || can_vertices_collide(fa1i, fb2i)
            || can_vertices_collide(fa2i, fb0i)
            || can_vertices_collide(fa2i, fb1i)
            || can_vertices_collide(fa2i, fb2i));
}

size_t LBVH::num_vertex_nodes() const
{
    return m_impl->vertex_bvh.nodes.size();
}

size_t LBVH::num_edge_nodes() const { return m_impl->edge_bvh.nodes.size(); }

size_t LBVH::num_face_nodes() const { return m_impl->face_bvh.nodes.size(); }

// ---------------------------------------------------------------------------
// Debug / validation.

void LBVH::vertex_nodes_to_host(
    ipc::LBVH::Nodes& nodes, ipc::LBVH::RightmostLeaves& rightmost_leaves) const
{
    to_host(m_impl->vertex_bvh, nodes, rightmost_leaves);
}

void LBVH::edge_nodes_to_host(
    ipc::LBVH::Nodes& nodes, ipc::LBVH::RightmostLeaves& rightmost_leaves) const
{
    to_host(m_impl->edge_bvh, nodes, rightmost_leaves);
}

void LBVH::face_nodes_to_host(
    ipc::LBVH::Nodes& nodes, ipc::LBVH::RightmostLeaves& rightmost_leaves) const
{
    to_host(m_impl->face_bvh, nodes, rightmost_leaves);
}

} // namespace ipc::cuda

#endif // IPC_TOOLKIT_WITH_CUDA
