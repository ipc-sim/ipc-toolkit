#include "lbvh.hpp"

#include <ipc/math/morton.hpp>
#include <ipc/utils/merge_thread_local.hpp>
#include <ipc/utils/profiler.hpp>

#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>

#if defined(__APPLE__)
// We utilize SIMD registers to compare 1 Node against 4 Queries simultaneously.
#include <simd/simd.h>
#endif

using namespace std::placeholders;

namespace ipc {

namespace {
    // Helper to safely convert double AABB to float AABB
    inline void assign_inflated_aabb(const AABB& box, LBVH::Node& node)
    {
        constexpr float inf = std::numeric_limits<float>::infinity();
        // Round Min down
        node.aabb_min[0] = std::nextafter(float(box.min[0]), -inf);
        node.aabb_min[1] = std::nextafter(float(box.min[1]), -inf);
        node.aabb_min[2] = box.min.size() == 3
            ? std::nextafter(float(box.min[2]), -inf)
            : 0.0f;
        // Round Max up
        node.aabb_max[0] = std::nextafter(float(box.max[0]), inf);
        node.aabb_max[1] = std::nextafter(float(box.max[1]), inf);
        node.aabb_max[2] =
            box.max.size() == 3 ? std::nextafter(float(box.max[2]), inf) : 0.0f;
    }
} // namespace

LBVH::LBVH() : BroadPhase() { }

LBVH::~LBVH() = default;

// Initialize defaults
LBVH::Node::Node() : primitive_id(INVALID_ID), is_inner_marker(0)
{
    static_assert(
        sizeof(Node) == 32,
        "LBVH::Node size must be 32 bytes to fit 2 Nodes in a cache line");
}

void LBVH::build(
    Eigen::ConstRef<Eigen::MatrixXi> edges,
    Eigen::ConstRef<Eigen::MatrixXi> faces)
{
    BroadPhase::build(edges, faces); // Build edge_boxes and face_boxes

    if (vertex_boxes.empty()) {
        return;
    }

    dim = vertex_boxes[0].min.size();
    assert(dim == 2 || dim == 3); // Only 2D and 3D supported

    mesh_aabb = { Eigen::Array3d::Zero(), Eigen::Array3d::Zero() };
    mesh_aabb.min = to_3D(vertex_boxes[0].min);
    mesh_aabb.max = to_3D(vertex_boxes[0].max);
    for (const auto& box : vertex_boxes) {
        mesh_aabb.min = mesh_aabb.min.min(to_3D(box.min));
        mesh_aabb.max = mesh_aabb.max.max(to_3D(box.max));
    }

    init_bvh(vertex_boxes, vertex_bvh);
    init_bvh(edge_boxes, edge_bvh);
    init_bvh(face_boxes, face_bvh);

    // Copy edge and face vertex ids for access during traversal
    edge_vertex_ids.resize(edges.rows());
    for (int i = 0; i < edges.rows(); i++) {
        edge_vertex_ids[i][0] = static_cast<index_t>(edges(i, 0));
        edge_vertex_ids[i][1] = static_cast<index_t>(edges(i, 1));
    }
    face_vertex_ids.resize(faces.rows());
    for (int i = 0; i < faces.rows(); i++) {
        face_vertex_ids[i][0] = static_cast<index_t>(faces(i, 0));
        face_vertex_ids[i][1] = static_cast<index_t>(faces(i, 1));
        face_vertex_ids[i][2] = static_cast<index_t>(faces(i, 2));
    }

    // Clear parent data to save memory.
    // These are redundant after building the BVHs.
    vertex_boxes.clear();
    vertex_boxes.shrink_to_fit();
    edge_boxes.clear();
    edge_boxes.shrink_to_fit();
    face_boxes.clear();
    face_boxes.shrink_to_fit();
}

namespace {

    int delta(
        const std::vector<LBVH::MortonCodeElement>& sorted_morton_codes,
        int i,
        uint64_t code_i,
        int j)
    {
        if (j < 0 || j > sorted_morton_codes.size() - 1) {
            return -1;
        }
        uint64_t code_j = sorted_morton_codes[j].morton_code;
        if (code_i == code_j) {
            // handle duplicate morton codes
            int element_idx_i = i; // sorted_morton_codes[i].elementIdx;
            int element_idx_j = j; // sorted_morton_codes[j].elementIdx;

            // add 32 for common prefix of code_i ^ code_j
#if defined(__GNUC__) || defined(__clang__)
            return 32 + __builtin_clz(element_idx_i ^ element_idx_j);
#elif defined(WIN32)
            return 32 + __lzcnt(element_idx_i ^ element_idx_j);
#endif
        }
#if defined(__GNUC__) || defined(__clang__)
        return __builtin_clzll(code_i ^ code_j);
#elif defined(WIN32)
        return __lzcnt64(code_i ^ code_j);
#endif
    }

    void determine_range(
        const std::vector<LBVH::MortonCodeElement>& sorted_morton_codes,
        int idx,
        int& lower,
        int& upper)
    {
        // determine direction of the range (+1 or -1)
        const uint64_t code = sorted_morton_codes[idx].morton_code;
        const int delta_l = delta(sorted_morton_codes, idx, code, idx - 1);
        const int delta_r = delta(sorted_morton_codes, idx, code, idx + 1);
        const int d = (delta_r >= delta_l) ? 1 : -1;

        // compute upper bound for the length of the range
        const int delta_min = std::min(delta_l, delta_r);
        int l_max = 2;
        while (delta(sorted_morton_codes, idx, code, idx + l_max * d)
               > delta_min) {
            l_max = l_max << 1;
        }

        // find the other end using binary search
        int l = 0;
        for (int t = l_max >> 1; t > 0; t >>= 1) {
            if (delta(sorted_morton_codes, idx, code, idx + (l + t) * d)
                > delta_min) {
                l += t;
            }
        }
        int jdx = idx + l * d;

        // ensure idx < jdx
        lower = std::min(idx, jdx);
        upper = std::max(idx, jdx);
    }

    int find_split(
        const std::vector<LBVH::MortonCodeElement>& sorted_morton_codes,
        int first,
        int last)
    {
        uint64_t first_code = sorted_morton_codes[first].morton_code;

        // Calculate the number of highest bits that are the same
        // for all objects, using the count-leading-zeros intrinsic.
        int common_prefix = delta(sorted_morton_codes, first, first_code, last);

        // Use binary search to find where the next bit differs.
        // Specifically, we are looking for the highest object that
        // shares more than common_prefix bits with the first one.
        int split = first; // initial guess
        int stride = last - first;
        do {
            stride = (stride + 1) >> 1;     // exponential decrease
            int new_split = split + stride; // proposed new position
            if (new_split < last) {
                int split_prefix =
                    delta(sorted_morton_codes, first, first_code, new_split);
                if (split_prefix > common_prefix) {
                    split = new_split; // accept proposal
                }
            }
        } while (stride > 1);

        return split;
    }
} // namespace

void LBVH::init_bvh(
    const std::vector<AABB>& boxes, std::vector<Node>& lbvh) const
{
    if (boxes.empty()) {
        return;
    }

    IPC_TOOLKIT_PROFILE_BLOCK("LBVH::init_bvh");

    std::vector<MortonCodeElement> morton_codes(boxes.size());
    {
        IPC_TOOLKIT_PROFILE_BLOCK("compute_morton_codes");
        const Eigen::Array3d mesh_width = mesh_aabb.max - mesh_aabb.min;
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, boxes.size()),
            [&](const tbb::blocked_range<size_t>& r) {
                for (size_t i = r.begin(); i < r.end(); i++) {
                    const auto& box = boxes[i];

                    const Eigen::Array3d center = 0.5 * (box.min + box.max);
                    const Eigen::Array3d mapped_center =
                        (center - mesh_aabb.min) / mesh_width;

                    if (dim == 2) {
                        morton_codes[i].morton_code =
                            morton_2D(mapped_center.x(), mapped_center.y());
                    } else {
                        morton_codes[i].morton_code = morton_3D(
                            mapped_center.x(), mapped_center.y(),
                            mapped_center.z());
                    }
                    morton_codes[i].box_id = i;
                }
            });
    }

    {
        IPC_TOOLKIT_PROFILE_BLOCK("sort_morton_codes");
        tbb::parallel_sort(
            morton_codes.begin(), morton_codes.end(),
            [](const MortonCodeElement& a, const MortonCodeElement& b) {
                return a.morton_code < b.morton_code;
            });
    }

    const int LEAF_OFFSET = int(boxes.size()) - 1;

    lbvh.resize(2 * boxes.size() - 1);
    std::vector<ConstructionInfo> construction_infos(2 * boxes.size() - 1);
    {
        IPC_TOOLKIT_PROFILE_BLOCK("build_hierarchy");
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, boxes.size()),
            [&](const tbb::blocked_range<size_t>& r) {
                for (size_t i = r.begin(); i < r.end(); i++) {

                    assert(i < boxes.size());
                    {
                        const auto& box = boxes[morton_codes[i].box_id];

                        Node leaf_node; // Create leaf node
                        assign_inflated_aabb(box, leaf_node);
                        leaf_node.primitive_id = morton_codes[i].box_id;
                        lbvh[LEAF_OFFSET + i] = leaf_node; // Store leaf
                    }

                    if (i < boxes.size() - 1) {
                        // Find out which range of objects the node corresponds
                        // to. (This is where the magic happens!)

                        int first, last;
                        determine_range(morton_codes, int(i), first, last);

                        // Determine where to split the range
                        int split = find_split(morton_codes, first, last);

                        // Select child_a
                        int child_a = -1;
                        if (split == first) {
                            // pointer to leaf node
                            child_a = LEAF_OFFSET + split;
                        } else {
                            child_a = split; // pointer to internal node
                        }

                        // Select child_b
                        int child_b = -1;
                        if (split + 1 == last) {
                            child_b =
                                LEAF_OFFSET + split + 1; // pointer to leaf node
                        } else {
                            child_b = split + 1; // pointer to internal node
                        }

                        // Record parent-child relationships
                        lbvh[i].left = child_a;
                        lbvh[i].right = child_b;
                        construction_infos[child_a].parent = int(i);
                        construction_infos[child_b].parent = int(i);
                    }

                    // node 0 is the root
                    if (i == 0) {
                        construction_infos[0].parent = 0;
                    }
                }
            });
    }

    {
        IPC_TOOLKIT_PROFILE_BLOCK("populate_boxes");
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, boxes.size()),
            [&](const tbb::blocked_range<size_t>& r) {
                for (size_t i = r.begin(); i < r.end(); i++) {
                    int node_idx = construction_infos[LEAF_OFFSET + i].parent;
                    while (true) {
                        auto& info = construction_infos[node_idx];
                        if (info.visitation_count++ == 0) {
                            // this is the first thread that arrived at this
                            // node -> finished
                            break;
                        }
                        // this is the second thread that arrived at this node,
                        // both children are computed -> compute aabb union and
                        // continue
                        assert(lbvh[node_idx].is_inner());
                        const Node& child_b = lbvh[lbvh[node_idx].right];
                        const Node& child_a = lbvh[lbvh[node_idx].left];
                        lbvh[node_idx].aabb_min =
                            child_a.aabb_min.min(child_b.aabb_min);
                        lbvh[node_idx].aabb_max =
                            child_a.aabb_max.max(child_b.aabb_max);

                        if (node_idx == 0) {
                            break; // root node
                        }
                        node_idx = info.parent;
                    }
                }
            });
    }
}

void LBVH::clear()
{
    BroadPhase::clear();
    // Clear BVH nodes
    vertex_bvh.clear();
    edge_bvh.clear();
    face_bvh.clear();

    // Clear vertex IDs
    edge_vertex_ids.clear();
    face_vertex_ids.clear();
}

namespace {
    template <typename Candidate, bool swap_order, bool triangular>
    inline void attempt_add_candidate(
        const LBVH::Node& query,
        const LBVH::Node& node,
        const std::function<bool(size_t, size_t)>& can_collide,
        std::vector<Candidate>& candidates)
    {
        int i = query.primitive_id, j = node.primitive_id;
        if constexpr (swap_order) {
            std::swap(i, j);
        }

        if constexpr (triangular) {
            if (i >= j) {
                return;
            }
        }

        if (!can_collide(i, j)) {
            return;
        }

        candidates.emplace_back(i, j);
    }

    template <typename Candidate, bool swap_order, bool triangular>
    void traverse_lbvh(
        const LBVH::Node& query,
        const std::vector<LBVH::Node>& lbvh,
        const std::function<bool(size_t, size_t)>& can_collide,
        std::vector<Candidate>& candidates)
    {
        // Use a fixed-size array as a stack to avoid dynamic allocations
        constexpr int MAX_STACK_SIZE = 64;
        int stack[MAX_STACK_SIZE];
        int stack_ptr = 0;
        stack[stack_ptr++] = LBVH::Node::INVALID_POINTER;

        int node_idx = 0; // root
        do {
            const LBVH::Node& node = lbvh[node_idx];

            // Check left and right are valid pointers
            assert(node.is_inner());

#if defined(__GNUC__) || defined(__clang__)
            // Prefetch child nodes to reduce cache misses
            __builtin_prefetch(&lbvh[node.left], 0, 1);
            __builtin_prefetch(&lbvh[node.right], 0, 1);
#endif

            const LBVH::Node& child_l = lbvh[node.left];
            const LBVH::Node& child_r = lbvh[node.right];
            const bool intersects_l = child_l.intersects(query);
            const bool intersects_r = child_r.intersects(query);

            // Query overlaps a leaf node => report collision.
            if (intersects_l && child_l.is_leaf()) {
                attempt_add_candidate<Candidate, swap_order, triangular>(
                    query, child_l, can_collide, candidates);
            }
            if (intersects_r && child_r.is_leaf()) {
                attempt_add_candidate<Candidate, swap_order, triangular>(
                    query, child_r, can_collide, candidates);
            }

            // Query overlaps an internal node => traverse.
            bool traverse_l = (intersects_l && !child_l.is_leaf());
            bool traverse_r = (intersects_r && !child_r.is_leaf());

            if (!traverse_l && !traverse_r) {
                assert(stack_ptr > 0);
                node_idx = stack[--stack_ptr];
            } else {
                node_idx = traverse_l ? node.left : node.right;
                if (traverse_l && traverse_r) {
                    // Postpone traversal of the right child
                    assert(stack_ptr < MAX_STACK_SIZE);
                    stack[stack_ptr++] = node.right;
                }
            }
        } while (node_idx != LBVH::Node::INVALID_POINTER); // Same as root
    }

#ifdef __APPLE__
    // SIMD Traversal
    // Traverses 4 queries simultaneously using SIMD.
    template <typename Candidate, bool swap_order, bool triangular>
    void traverse_lbvh_simd(
        const LBVH::Node* queries,
        const size_t n_queries,
        const std::vector<LBVH::Node>& lbvh,
        const std::function<bool(size_t, size_t)>& can_collide,
        std::vector<Candidate>& candidates)
    {
        assert(n_queries >= 1 && n_queries <= 4);
        // Load 4 queries into single registers (Structure of Arrays)
        auto make_simd = [&](auto F) -> simd_float4 {
            return simd_float4 {
                F(0),
                n_queries > 1 ? F(1) : 0.0f,
                n_queries > 2 ? F(2) : 0.0f,
                n_queries > 3 ? F(3) : 0.0f,
            };
        };

        const simd_float4 q_min_x =
            make_simd([&](int k) { return queries[k].aabb_min.x(); });
        const simd_float4 q_min_y =
            make_simd([&](int k) { return queries[k].aabb_min.y(); });
        const simd_float4 q_min_z =
            make_simd([&](int k) { return queries[k].aabb_min.z(); });
        const simd_float4 q_max_x =
            make_simd([&](int k) { return queries[k].aabb_max.x(); });
        const simd_float4 q_max_y =
            make_simd([&](int k) { return queries[k].aabb_max.y(); });
        const simd_float4 q_max_z =
            make_simd([&](int k) { return queries[k].aabb_max.z(); });

        // Use a fixed-size array as a stack to avoid dynamic allocations
        constexpr int MAX_STACK_SIZE = 64;
        int stack[MAX_STACK_SIZE];
        int stack_ptr = 0;
        stack[stack_ptr++] = LBVH::Node::INVALID_POINTER;

        int node_idx = 0; // root
        do {
            const LBVH::Node& node = lbvh[node_idx];

            // Check left and right are valid pointers
            assert(node.is_inner());

#if defined(__GNUC__) || defined(__clang__)
            // Prefetch child nodes to reduce cache misses
            __builtin_prefetch(&lbvh[node.left], 0, 1);
            __builtin_prefetch(&lbvh[node.right], 0, 1);
#endif

            const LBVH::Node& child_l = lbvh[node.left];
            const LBVH::Node& child_r = lbvh[node.right];

            // 1. Intersect 4 queries at once
            // (child_l.min <= query.max) && (query.min <= child_l.max)
            const simd_int4 intersects_l = (child_l.aabb_min.x() <= q_max_x)
                & (child_l.aabb_min.y() <= q_max_y)
                & (child_l.aabb_min.z() <= q_max_z)
                & (q_min_x <= child_l.aabb_max.x())
                & (q_min_y <= child_l.aabb_max.y())
                & (q_min_z <= child_l.aabb_max.z());

            // 2. Intersect 4 queries at once
            // (child_r.min <= query.max) && (query.min <= child_r.max)
            const simd_int4 intersects_r = (child_r.aabb_min.x() <= q_max_x)
                & (child_r.aabb_min.y() <= q_max_y)
                & (child_r.aabb_min.z() <= q_max_z)
                & (q_min_x <= child_r.aabb_max.x())
                & (q_min_y <= child_r.aabb_max.y())
                & (q_min_z <= child_r.aabb_max.z());

            const bool any_intersects_l = simd_any(intersects_l);
            const bool any_intersects_r = simd_any(intersects_r);

            // Query overlaps a leaf node => report collision
            if (any_intersects_l && child_l.is_leaf()) {
                for (int k = 0; k < n_queries; ++k) {
                    if (intersects_l[k]) {
                        attempt_add_candidate<
                            Candidate, swap_order, triangular>(
                            queries[k], child_l, can_collide, candidates);
                    }
                }
            }
            if (any_intersects_r && child_r.is_leaf()) {
                for (int k = 0; k < n_queries; ++k) {
                    if (intersects_r[k]) {
                        attempt_add_candidate<
                            Candidate, swap_order, triangular>(
                            queries[k], child_r, can_collide, candidates);
                    }
                }
            }

            // Query overlaps an internal node => traverse.
            bool traverse_l = (any_intersects_l && !child_l.is_leaf());
            bool traverse_r = (any_intersects_r && !child_r.is_leaf());

            if (!traverse_l && !traverse_r) {
                assert(stack_ptr > 0);
                node_idx = stack[--stack_ptr];
            } else {
                node_idx = traverse_l ? node.left : node.right;
                if (traverse_l && traverse_r) {
                    // Postpone traversal of the right child
                    assert(stack_ptr < MAX_STACK_SIZE);
                    stack[stack_ptr++] = node.right;
                }
            }
        } while (node_idx != LBVH::Node::INVALID_POINTER); // Same as root
    }
#endif

    template <
        typename Candidate,
        bool swap_order,
        bool triangular,
        bool use_simd = true>
    void independent_traversal(
        const std::vector<LBVH::Node>& source,
        const std::vector<LBVH::Node>& target,
        const std::function<bool(size_t, size_t)>& can_collide,
        tbb::enumerable_thread_specific<std::vector<Candidate>>& storage)
    {
#ifdef __APPLE__ // Only support SIMD on Apple platforms for now
        constexpr size_t SIMD_SIZE = use_simd ? 4 : 1;
        constexpr size_t GRAIN_SIZE = use_simd ? 64 : 1;
#else
        constexpr size_t SIMD_SIZE = 1;
        constexpr size_t GRAIN_SIZE = 1;
#endif

        // Calculate the offset to the first leaf node in the source BVH.
        const size_t SOURCE_LEAF_OFFSET = source.size() / 2;
        const size_t N_SOURCE_LEAVES = SOURCE_LEAF_OFFSET + 1;

        tbb::parallel_for(
            tbb::blocked_range<size_t>(size_t(0), N_SOURCE_LEAVES, GRAIN_SIZE),
            [&](const tbb::blocked_range<size_t>& r) {
                auto& local_candidates = storage.local();
                for (size_t i = r.begin(); i < r.end(); i += SIMD_SIZE) {
#ifdef __APPLE__
                    if constexpr (use_simd) {
                        traverse_lbvh_simd<Candidate, swap_order, triangular>(
                            &source[SOURCE_LEAF_OFFSET + i],
                            std::min(SIMD_SIZE, r.end() - i), target,
                            can_collide, local_candidates);
                    } else {
#endif
                        traverse_lbvh<Candidate, swap_order, triangular>(
                            source[SOURCE_LEAF_OFFSET + i], target, can_collide,
                            local_candidates);
#ifdef __APPLE__
                    }
#endif
                }
            });
    }
} // namespace

template <typename Candidate, bool swap_order, bool triangular>
void LBVH::detect_candidates(
    const std::vector<Node>& source,
    const std::vector<Node>& target,
    const std::function<bool(size_t, size_t)>& can_collide,
    std::vector<Candidate>& candidates)
{
    if (source.empty() || target.empty()) {
        return;
    }

    tbb::enumerable_thread_specific<std::vector<Candidate>> storage;

    independent_traversal<Candidate, swap_order, triangular>(
        source, target, can_collide, storage);

    merge_thread_local_vectors(storage, candidates);
}

void LBVH::detect_vertex_vertex_candidates(
    std::vector<VertexVertexCandidate>& candidates) const
{
    if (!has_vertices()) {
        return;
    }

    IPC_TOOLKIT_PROFILE_BLOCK("LBVH::detect_vertex_vertex_candidates");

    detect_candidates(vertex_bvh, can_vertices_collide, candidates);
}

void LBVH::detect_edge_vertex_candidates(
    std::vector<EdgeVertexCandidate>& candidates) const
{
    if (!has_edges() || !has_vertices()) {
        return;
    }

    IPC_TOOLKIT_PROFILE_BLOCK("LBVH::detect_edge_vertex_candidates");

    // In 2D and for codimensional edge-vertex collisions, there are more
    // vertices than edges, so we want to iterate over the edges.
    detect_candidates(
        edge_bvh, vertex_bvh,
        std::bind(&LBVH::can_edge_vertex_collide, this, _1, _2), candidates);
}

void LBVH::detect_edge_edge_candidates(
    std::vector<EdgeEdgeCandidate>& candidates) const
{
    if (!has_edges()) {
        return;
    }

    IPC_TOOLKIT_PROFILE_BLOCK("LBVH::detect_edge_edge_candidates");

    detect_candidates(
        edge_bvh, std::bind(&LBVH::can_edges_collide, this, _1, _2),
        candidates);
}

void LBVH::detect_face_vertex_candidates(
    std::vector<FaceVertexCandidate>& candidates) const
{
    if (!has_faces() || !has_vertices()) {
        return;
    }

    IPC_TOOLKIT_PROFILE_BLOCK("LBVH::detect_face_vertex_candidates");

    // The ratio vertices:faces is 1:2, so we want to iterate over the vertices.
    detect_candidates<FaceVertexCandidate, /*swap_order=*/true>(
        vertex_bvh, face_bvh,
        std::bind(&LBVH::can_face_vertex_collide, this, _1, _2), candidates);
}

void LBVH::detect_edge_face_candidates(
    std::vector<EdgeFaceCandidate>& candidates) const
{
    if (!has_edges() || !has_faces()) {
        return;
    }

    IPC_TOOLKIT_PROFILE_BLOCK("LBVH::detect_edge_face_candidates");

    // The ratio edges:faces is 3:2, so we want to iterate over the faces.
    detect_candidates<EdgeFaceCandidate, /*swap_order=*/true>(
        face_bvh, edge_bvh,
        std::bind(&LBVH::can_edge_face_collide, this, _1, _2), candidates);
}

void LBVH::detect_face_face_candidates(
    std::vector<FaceFaceCandidate>& candidates) const
{
    if (!has_faces()) {
        return;
    }

    IPC_TOOLKIT_PROFILE_BLOCK("LBVH::detect_face_face_candidates");
    detect_candidates(
        face_bvh, std::bind(&LBVH::can_faces_collide, this, _1, _2),
        candidates);
}

// ============================================================================

bool LBVH::can_edge_vertex_collide(size_t ei, size_t vi) const
{
    assert(ei < edge_vertex_ids.size());
    const auto& [e0i, e1i] = edge_vertex_ids[ei];

    return vi != e0i && vi != e1i
        && (can_vertices_collide(vi, e0i) || can_vertices_collide(vi, e1i));
}

bool LBVH::can_edges_collide(size_t eai, size_t ebi) const
{
    assert(eai < edge_vertex_ids.size());
    const auto& [ea0i, ea1i] = edge_vertex_ids[eai];
    assert(ebi < edge_vertex_ids.size());
    const auto& [eb0i, eb1i] = edge_vertex_ids[ebi];

    const bool share_endpoint =
        ea0i == eb0i || ea0i == eb1i || ea1i == eb0i || ea1i == eb1i;

    return !share_endpoint
        && (can_vertices_collide(ea0i, eb0i) || can_vertices_collide(ea0i, eb1i)
            || can_vertices_collide(ea1i, eb0i)
            || can_vertices_collide(ea1i, eb1i));
}

bool LBVH::can_face_vertex_collide(size_t fi, size_t vi) const
{
    assert(fi < face_vertex_ids.size());
    const auto& [f0i, f1i, f2i] = face_vertex_ids[fi];

    return vi != f0i && vi != f1i && vi != f2i
        && (can_vertices_collide(vi, f0i) || can_vertices_collide(vi, f1i)
            || can_vertices_collide(vi, f2i));
}

bool LBVH::can_edge_face_collide(size_t ei, size_t fi) const
{
    assert(ei < edge_vertex_ids.size());
    const auto& [e0i, e1i] = edge_vertex_ids[ei];
    assert(fi < face_vertex_ids.size());
    const auto& [f0i, f1i, f2i] = face_vertex_ids[fi];

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
    assert(fai < face_vertex_ids.size());
    const auto& [fa0i, fa1i, fa2i] = face_vertex_ids[fai];
    assert(fbi < face_vertex_ids.size());
    const auto& [fb0i, fb1i, fb2i] = face_vertex_ids[fbi];

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

} // namespace ipc