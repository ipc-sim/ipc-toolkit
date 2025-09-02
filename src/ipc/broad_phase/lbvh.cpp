#include "lbvh.hpp"

#include <ipc/utils/merge_thread_local.hpp>
#include <ipc/utils/profiler.hpp>

#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>

using namespace std::placeholders;

namespace ipc {

namespace {
    // Expands a 21-bit integer into 63 bits by inserting 2 zeros after each
    // bit.
    uint64_t expand_bits(uint64_t v)
    {
        v = (v | v << 32) & 0x1F00000000FFFF;
        v = (v | v << 16) & 0x1F0000FF0000FF;
        v = (v | v << 8) & 0x100F00F00F00F00F;
        v = (v | v << 4) & 0x10C30C30C30C30C3;
        v = (v | v << 2) & 0x1249249249249249;
        return v;
    }

    // Calculates a 63-bit Morton code for the given 3D point located within the
    // unit cube [0,1].
    uint64_t morton_3D(double x, double y, double z)
    {
        constexpr double scale = 1 << 21;
        x = std::clamp(x * scale, 0.0, scale - 1);
        y = std::clamp(y * scale, 0.0, scale - 1);
        z = std::clamp(z * scale, 0.0, scale - 1);
        uint64_t xx = expand_bits(uint64_t(x));
        uint64_t yy = expand_bits(uint64_t(y));
        uint64_t zz = expand_bits(uint64_t(z));
        return (xx << 2) | (yy << 1) | zz;
    }
} // namespace

LBVH::LBVH() : BroadPhase() { }

LBVH::~LBVH() = default;

void LBVH::build(
    Eigen::ConstRef<Eigen::MatrixXi> edges,
    Eigen::ConstRef<Eigen::MatrixXi> faces)
{
    BroadPhase::build(edges, faces); // Build edge_boxes and face_boxes

    assert(vertex_boxes.size() > 0);
    mesh_aabb = { vertex_boxes[0].min, vertex_boxes[0].max };
    for (const auto& box : vertex_boxes) {
        mesh_aabb.min = mesh_aabb.min.min(box.min);
        mesh_aabb.max = mesh_aabb.max.max(box.max);
    }

    init_bvh(vertex_boxes, vertex_bvh);
    init_bvh(edge_boxes, edge_bvh);
    init_bvh(face_boxes, face_bvh);
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
    if (boxes.size() == 0) {
        return;
    }

    IPC_TOOLKIT_PROFILE_BLOCK("LBVH::init_bvh");

    std::vector<MortonCodeElement> morton_codes(boxes.size());
    {
        IPC_TOOLKIT_PROFILE_BLOCK("compute_morton_codes");
        const ArrayMax3d mesh_width = mesh_aabb.max - mesh_aabb.min;
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, boxes.size()),
            [&](const tbb::blocked_range<size_t>& r) {
                for (size_t i = r.begin(); i < r.end(); i++) {
                    const auto& box = boxes[i];

                    const ArrayMax3d center = 0.5 * (box.min + box.max);
                    const ArrayMax3d mapped_center =
                        (center - mesh_aabb.min) / mesh_width;

                    morton_codes[i].morton_code = morton_3D(
                        mapped_center.x(), mapped_center.y(),
                        mapped_center.z());
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
                        lbvh[LEAF_OFFSET + i].left = Node::INVALID_POINTER;
                        lbvh[LEAF_OFFSET + i].right = Node::INVALID_POINTER;
                        lbvh[LEAF_OFFSET + i].aabb_min = box.min;
                        lbvh[LEAF_OFFSET + i].aabb_max = box.max;
                        lbvh[LEAF_OFFSET + i].vertex_ids = box.vertex_ids;
                        lbvh[LEAF_OFFSET + i].primitive_id =
                            morton_codes[i].box_id;
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
                        if constexpr (Node::ABSOLUTE_POINTERS) {
                            lbvh[i].left = child_a;
                            lbvh[i].right = child_b;
                        } else {
                            lbvh[i].left = child_a - int(i);
                            lbvh[i].right = child_b - int(i);
                        }
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
                        const Node& child_b =
                            lbvh[Node::POINTER(node_idx, lbvh[node_idx].right)];
                        const Node& child_a =
                            lbvh[Node::POINTER(node_idx, lbvh[node_idx].left)];

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
    vertex_bvh.clear();
    edge_bvh.clear();
    face_bvh.clear();
}

namespace {
    inline bool share_a_vertex(
        const std::array<index_t, 3>& a, const std::array<index_t, 3>& b)
    {
        return a[0] == b[0] || a[0] == b[1] || a[0] == b[2]
            || (a[1] >= 0 && (a[1] == b[0] || a[1] == b[1] || a[1] == b[2]))
            || (a[2] >= 0 && (a[2] == b[0] || a[2] == b[1] || a[2] == b[2]));
    }

    template <typename Candidate, bool swap_order, bool triangular>
    inline void attempt_add_candidate(
        const LBVH::Node& query,
        const LBVH::Node& node,
        const std::function<bool(size_t, size_t)>& can_collide,
        std::vector<Candidate>& candidates)
    {
        if (share_a_vertex(query.vertex_ids, node.vertex_ids)) {
            return;
        }

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

            // Prefetch child nodes to reduce cache misses
            __builtin_prefetch(&lbvh[node.left], 0, 1);
            __builtin_prefetch(&lbvh[node.right], 0, 1);

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
                node_idx = stack[--stack_ptr];
            } else {
                node_idx = traverse_l ? node.left : node.right;
                if (traverse_l && traverse_r) {
                    // Postpone traversal of the right child
                    stack[stack_ptr++] = node.right;
                }
            }
        } while (node_idx != LBVH::Node::INVALID_POINTER); // Same as root
    }

    // Traverses the target BVH independently for each leaf node in the source
    // BVH. For each source leaf, performs a stack-based traversal of the target
    // BVH, collecting candidate pairs that pass the can_collide predicate.
    // Results are stored in thread-local storage for later merging.
    template <typename Candidate, bool swap_order, bool triangular>
    void independent_traversal(
        const std::vector<LBVH::Node>& source,
        const std::vector<LBVH::Node>& target,
        const std::function<bool(size_t, size_t)>& can_collide,
        tbb::enumerable_thread_specific<std::vector<Candidate>>& storage)
    {
        // Calculate the offset to the first leaf node in the source BVH.
        const size_t SOURCE_LEAF_OFFSET = source.size() / 2;
        const size_t N_SOURCE_LEAVES = SOURCE_LEAF_OFFSET + 1;

        tbb::parallel_for(
            tbb::blocked_range<size_t>(size_t(0), N_SOURCE_LEAVES),
            [&](const tbb::blocked_range<size_t>& r) {
                auto& local_candidates = storage.local();

                for (size_t i = r.begin(); i < r.end(); ++i) {
                    const auto& query_node = source[SOURCE_LEAF_OFFSET + i];
                    traverse_lbvh<Candidate, swap_order, triangular>(
                        query_node, target, can_collide, local_candidates);
                }
            });
    }

    // Parallel traversal of two BVHs using TBB task_group.
    // Recursively explores all pairs of nodes whose bounding boxes intersect,
    // adding candidate pairs when both nodes are leaves.
    template <typename Candidate, bool swap_order, bool triangular>
    void traverse_lbvh(
        const std::vector<LBVH::Node>& source,
        const std::vector<LBVH::Node>& target,
        int source_idx,
        int target_idx,
        const std::function<bool(size_t, size_t)>& can_collide,
        tbb::task_group& g,
        tbb::enumerable_thread_specific<std::vector<Candidate>>& storage)
    {
        // 1. Check for bounding box intersection
        if (!source[source_idx].intersects(target[target_idx])) {
            return;
        }

        // 2. Handle leaf nodes (base case)
        if (source[source_idx].is_leaf() && target[target_idx].is_leaf()) {
            attempt_add_candidate<Candidate, swap_order, triangular>(
                source[source_idx], target[target_idx], can_collide,
                storage.local());
            return;
        }

        // 3. Handle mixed or internal nodes (parallel recursion)

        // TBB's task_group provides an easy way to offload tasks.
        const auto dispatch = [&](int source_idx, int target_idx) {
            g.run([&, source_idx, target_idx] {
                traverse_lbvh<Candidate, swap_order, triangular>(
                    source, target, source_idx, target_idx, can_collide, g,
                    storage);
            });
        };

        if (source[source_idx].is_leaf()) {
            dispatch(source_idx, target[target_idx].left);
            dispatch(source_idx, target[target_idx].right);
        } else if (target[target_idx].is_leaf()) {
            dispatch(source[source_idx].left, target_idx);
            dispatch(source[source_idx].right, target_idx);
        } else {
            // Both internal nodes, test all four combinations.
            dispatch(source[source_idx].left, target[target_idx].left);
            dispatch(source[source_idx].left, target[target_idx].right);
            dispatch(source[source_idx].right, target[target_idx].left);
            dispatch(source[source_idx].right, target[target_idx].right);
        }
    }

    template <typename Candidate, bool swap_order, bool triangular>
    void simultaneous_traversal(
        const std::vector<LBVH::Node>& source,
        const std::vector<LBVH::Node>& target,
        const std::function<bool(size_t, size_t)>& can_collide,
        tbb::enumerable_thread_specific<std::vector<Candidate>>& storage)
    {
        tbb::task_group g; // TBB task group to manage parallel work

        g.run_and_wait([&] {
            traverse_lbvh<Candidate, swap_order, triangular>(
                source, target, /*source_root*/ 0, /*target_root*/ 0,
                can_collide, g, storage);
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
    tbb::enumerable_thread_specific<std::vector<Candidate>> storage;

    independent_traversal<Candidate, swap_order, triangular>(
        source, target, can_collide, storage);
    // simultaneous_traversal<Candidate, swap_order, triangular>(
    //     source, target, can_collide, storage);

    merge_thread_local_vectors(storage, candidates);
}

void LBVH::detect_vertex_vertex_candidates(
    std::vector<VertexVertexCandidate>& candidates) const
{
    if (vertex_boxes.size() == 0) {
        return;
    }

    IPC_TOOLKIT_PROFILE_BLOCK("LBVH::detect_vertex_vertex_candidates");

    detect_candidates(vertex_bvh, can_vertices_collide, candidates);
}

void LBVH::detect_edge_vertex_candidates(
    std::vector<EdgeVertexCandidate>& candidates) const
{
    if (edge_boxes.size() == 0 || vertex_boxes.size() == 0) {
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
    if (edge_boxes.size() == 0) {
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
    if (face_boxes.size() == 0 || vertex_boxes.size() == 0) {
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
    if (edge_boxes.size() == 0 || face_boxes.size() == 0) {
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
    if (face_boxes.size() == 0) {
        return;
    }

    IPC_TOOLKIT_PROFILE_BLOCK("LBVH::detect_face_face_candidates");
    detect_candidates(
        face_bvh, std::bind(&LBVH::can_faces_collide, this, _1, _2),
        candidates);
}
} // namespace ipc