#include "lbvh.hpp"

#include <ipc/math/morton.hpp>
#include <ipc/utils/merge_thread_local.hpp>
#include <ipc/utils/profiler.hpp>

#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>

#ifdef IPC_TOOLKIT_WITH_SIMD
// We utilize SIMD registers to compare one node against multiple queries
// simultaneously, with the number of queries determined by
// xs::batch<float>::size.
#include <xsimd/xsimd.hpp>
namespace xs = xsimd;
#endif

#include <array>
#include <atomic>

using namespace std::placeholders;

namespace ipc {

namespace {
    // Helper to safely convert double AABB to float AABB
    inline void assign_inflated_aabb(const AABB& box, LBVH::Node& node)
    {
        // Round Min down
        node.aabb_min = box.min.unaryExpr([](double val) {
            return std::nextafter(
                float(val), -std::numeric_limits<float>::infinity());
        });
        // Round Max up
        node.aabb_max = box.max.unaryExpr([](double val) {
            return std::nextafter(
                float(val), std::numeric_limits<float>::infinity());
        });
    }
} // namespace

LBVH::LBVH() : BroadPhase()
{
    static_assert(
        sizeof(LBVH::Node) == 32,
        "LBVH::Node size must be 32 bytes to fit 2 Nodes in a cache line");
}

LBVH::~LBVH() = default;

void LBVH::build(
    Eigen::ConstRef<Eigen::MatrixXi> edges,
    Eigen::ConstRef<Eigen::MatrixXi> faces)
{
    IPC_TOOLKIT_PROFILE_BLOCK("LBVH::build");

    BroadPhase::build(edges, faces); // Build edge_boxes and face_boxes

    if (vertex_boxes.empty()) {
        return;
    }

    assert(dim == 2 || dim == 3);

    compute_mesh_aabb(mesh_aabb.min, mesh_aabb.max);

    init_bvh(vertex_boxes, vertex_bvh, vertex_rightmost_leaves);
    init_bvh(edge_boxes, edge_bvh, edge_rightmost_leaves);
    init_bvh(face_boxes, face_bvh, face_rightmost_leaves);

    // Copy edge and face vertex ids for access during traversal
    {
        IPC_TOOLKIT_PROFILE_BLOCK("copy_vertex_ids");

        edge_vertex_ids.resize(edges.rows());
        tbb::parallel_for(size_t(0), size_t(edges.rows()), [&](size_t i) {
            edge_vertex_ids[i][0] = static_cast<index_t>(edges(i, 0));
            edge_vertex_ids[i][1] = static_cast<index_t>(edges(i, 1));
        });

        face_vertex_ids.resize(faces.rows());
        tbb::parallel_for(size_t(0), size_t(faces.rows()), [&](size_t i) {
            face_vertex_ids[i][0] = static_cast<index_t>(faces(i, 0));
            face_vertex_ids[i][1] = static_cast<index_t>(faces(i, 1));
            face_vertex_ids[i][2] = static_cast<index_t>(faces(i, 2));
        });
    }

    // Clear parent data to save memory.
    // These are redundant after building the BVHs.
    vertex_boxes.clear();
    edge_boxes.clear();
    face_boxes.clear();
}

namespace {
    /// Returns the number of common leading bits (CLZ of XOR) between sorted
    /// Morton codes at positions i and j. code_i is the Morton code at position
    /// i, passed explicitly to avoid a redundant lookup. Returns -1 when j is
    /// out of bounds.  Duplicate codes fall back to CLZ of the index XOR
    /// (offset by 32 so it sorts after any code-level difference).
    int delta(
        const LBVH::MortonCodeElements& sorted_morton_codes,
        int i,
        uint64_t code_i,
        int j)
    {
        if (j < 0 || j >= sorted_morton_codes.size()) {
            return -1;
        }
        uint64_t code_j = sorted_morton_codes[j].morton_code;
        if (code_i == code_j) {
            // handle duplicate morton codes
            int element_idx_i = i;
            int element_idx_j = j;

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
} // namespace

void LBVH::init_bvh(
    const AABBs& boxes, Nodes& lbvh, RightmostLeaves& rightmost_leaves) const
{
    if (boxes.empty()) {
        return;
    }

    IPC_TOOLKIT_PROFILE_BLOCK("LBVH::init_bvh");

    if (lbvh.size() != 2 * boxes.size() - 1) {
        IPC_TOOLKIT_PROFILE_BLOCK("resize_bvh");
        lbvh.resize(2 * boxes.size() - 1);
    }

    LBVH::MortonCodeElements morton_codes(boxes.size());
    {
        IPC_TOOLKIT_PROFILE_BLOCK("compute_morton_codes");

        const Eigen::Array3d mesh_width = mesh_aabb.max - mesh_aabb.min;
        tbb::parallel_for(size_t(0), boxes.size(), [&](size_t i) {
            const auto& box = boxes[i];

            const Eigen::Array3d center = 0.5 * (box.min + box.max);
            const Eigen::Array3d mapped_center =
                (center - mesh_aabb.min) / mesh_width;

            if (dim == 2) {
                morton_codes[i].morton_code =
                    morton_2D(mapped_center.x(), mapped_center.y());
            } else {
                morton_codes[i].morton_code = morton_3D(
                    mapped_center.x(), mapped_center.y(), mapped_center.z());
            }
            morton_codes[i].box_id = i;
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

    assert(boxes.size() <= std::numeric_limits<int>::max());
    const int N_LEAVES = int(boxes.size());
    const int LEAF_OFFSET = N_LEAVES - 1;

    if (rightmost_leaves.size() != lbvh.size()) {
        rightmost_leaves.resize(lbvh.size());
    }

    LBVH::ConstructionInfos construction_infos(lbvh.size());
    {
        IPC_TOOLKIT_PROFILE_BLOCK("init_visitation_counts");
        tbb::parallel_for(size_t(0), lbvh.size(), [&](size_t i) {
            construction_infos[i].visitation_count.store(
                0, std::memory_order_relaxed);
        });
    }

    // Apetrei 2014: single bottom-up pass that simultaneously builds the
    // hierarchy and computes bounding boxes. Each leaf thread walks toward the
    // root, choosing its parent in O(1) by comparing the CLZ-delta values at
    // the two ends of its current key range.
    //
    // In this layout internal node j always splits between sorted keys j and
    // j+1. The root is NOT necessarily at index 0, so after construction we
    // swap the root into position 0 to match the traversal code's expectation.
    std::atomic<int> root_idx(-1);
    {
        IPC_TOOLKIT_PROFILE_BLOCK("build_hierarchy_and_boxes");
        tbb::parallel_for(0, N_LEAVES, [&](int i) {
            // --- Initialize leaf node ---
            {
                const auto& box = boxes[morton_codes[i].box_id];

                Node leaf_node; // Create leaf node
                assign_inflated_aabb(box, leaf_node);
                leaf_node.primitive_id = morton_codes[i].box_id;
                leaf_node.is_inner_marker = 0;
                lbvh[LEAF_OFFSET + i] = leaf_node; // Store leaf
                // A leaf's rightmost leaf is itself
                rightmost_leaves[LEAF_OFFSET + i] = i;
            }

            // --- Bottom-up walk (Apetrei 2014, Fig. 2) ---
            // Invariant: the current subtree covers the sorted-key range
            // [left_key, right_key].
            int left_key = i;
            int right_key = i;
            int current_node = LEAF_OFFSET + i;

            while (true) {
                // Choose parent. Candidates are internal node right_key
                // (current becomes its left / childA) or internal node
                // left_key-1 (current becomes its right / childB). Our delta()
                // returns CLZ (higher = more-similar = finer split), so the
                // CLOSER ancestor has the LARGER delta — hence ">".
                //
                // Boundary rules:
                //   left_key == 0        → must be childA (no node -1)
                //   right_key == n-1     → must be childB (no node n-1)
                const bool is_child_a = (left_key == 0)
                    || (right_key != N_LEAVES - 1
                        && delta(
                               morton_codes, right_key,
                               morton_codes[right_key].morton_code,
                               right_key + 1)
                            > delta(
                                morton_codes, left_key - 1,
                                morton_codes[left_key - 1].morton_code,
                                left_key));
                const int parent = is_child_a ? right_key : left_key - 1;

                auto& info = construction_infos[parent];

                // Write the child pointer on the parent node.
                // childA writes .left; childB writes .right.
                if (is_child_a) {
                    lbvh[parent].left = current_node;
                    info.left_range = left_key;
                } else {
                    lbvh[parent].right = current_node;
                    info.right_range = right_key;
                }

                // Atomic arrival gate: the first thread to reach this parent
                // stops; the second thread proceeds (it now knows both children
                // are complete).

                if (info.visitation_count++ == 0) {
                    // this is the first thread that arrived at this
                    // node -> finished
                    break;
                }
                // this is the second thread that arrived at this node,
                // both children are computed -> compute aabb union and
                // continue
                assert(lbvh[parent].is_inner());
                const Node& child_a = lbvh[lbvh[parent].left];
                const Node& child_b = lbvh[lbvh[parent].right];
                lbvh[parent].aabb_min = child_a.aabb_min.min(child_b.aabb_min);
                lbvh[parent].aabb_max = child_a.aabb_max.max(child_b.aabb_max);

                // Compute rightmost leaf: max of children's rightmost
                rightmost_leaves[parent] = std::max(
                    rightmost_leaves[lbvh[parent].left],
                    rightmost_leaves[lbvh[parent].right]);

                // Reconstruct the full key range for the parent.
                left_key = construction_infos[parent].left_range;
                right_key = construction_infos[parent].right_range;
                current_node = parent;

                if (left_key == 0 && right_key == N_LEAVES - 1) {
                    // only one thread should reach the root
                    int expected = -1;
                    [[maybe_unused]] bool set =
                        root_idx.compare_exchange_strong(
                            expected, current_node);
                    assert(set);
                    break; // root AABB is complete
                }
            }
        });
    }

    // --- Move the root to index 0 so traversal can start there. ---
    // In the Apetrei layout the root's index equals the global split position,
    // which is generally != 0.  We swap the root node into position 0 and patch
    // up the single affected child pointer.
    //
    // Key invariant (Apetrei): node 0's subtree always has left_key=0, so it is
    // only ever written as a LEFT child — meaning no internal node ever has
    // right==0.  Therefore swapping node 0 cannot create a spurious
    // is_inner_marker==0 (which would look like a leaf).
    const int root = root_idx.load();
    if (root > 0) {
        IPC_TOOLKIT_PROFILE_BLOCK("swap_root_to_zero");
        std::swap(lbvh[0], lbvh[root]);
        std::swap(rightmost_leaves[0], rightmost_leaves[root]);

        // The root (now at 0) is never any node's child, so no pointer
        // references R that needs rewriting to 0. The only pointers that
        // referenced 0 (the old node-0) must be rewritten to R.  And since old
        // node-0 was only ever a LEFT child (see invariant above), we only need
        // to patch .left pointers.
        tbb::parallel_for(size_t(0), lbvh.size(), [&](size_t i) {
            if (lbvh[i].is_inner() && lbvh[i].left == 0) {
                lbvh[i].left = root;
            }
        });
    }
}

void LBVH::clear()
{
    BroadPhase::clear();
    // Clear BVH nodes
    vertex_bvh.clear();
    vertex_rightmost_leaves.clear();
    edge_bvh.clear();
    edge_rightmost_leaves.clear();
    face_bvh.clear();
    face_rightmost_leaves.clear();

    // Clear vertex IDs
    edge_vertex_ids.clear();
    face_vertex_ids.clear();
}

namespace {
    template <typename Candidate, bool swap_order>
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

        if (!can_collide(i, j)) {
            return;
        }

        candidates.emplace_back(i, j);
    }

    template <typename Candidate, bool swap_order, bool triangular>
    void traverse_lbvh(
        const LBVH::Node& query,
        const size_t query_leaf_idx,
        const LBVH::Nodes& lbvh,
        const LBVH::RightmostLeaves& rightmost_leaves,
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

            if (lbvh.size() == 1) {     // Single node case (only root)
                assert(node.is_leaf()); // Only one node, so it must be a leaf
                if constexpr (triangular) {
                    break; // No self-collision if only one node
                }
                if (node.intersects(query)) {
                    attempt_add_candidate<Candidate, swap_order>(
                        query, node, can_collide, candidates);
                }
                break;
            }

            // Check left and right are valid pointers
            assert(node.is_inner());

#if defined(__GNUC__) || defined(__clang__)
            // Prefetch child nodes to reduce cache misses
            __builtin_prefetch(&lbvh[node.left], 0, 1);
            __builtin_prefetch(&lbvh[node.right], 0, 1);
#endif

            const LBVH::Node& child_l = lbvh[node.left];
            const LBVH::Node& child_r = lbvh[node.right];
            bool intersects_l = child_l.intersects(query);
            bool intersects_r = child_r.intersects(query);

            // Ignore overlap if the subtree is fully on the
            // left-hand side of the query (triangular traversal only).
            if constexpr (triangular) {
                if (intersects_l
                    && rightmost_leaves[node.left] <= query_leaf_idx) {
                    intersects_l = false;
                }
                if (intersects_r
                    && rightmost_leaves[node.right] <= query_leaf_idx) {
                    intersects_r = false;
                }
            }

            // Query overlaps a leaf node => report collision.
            if (intersects_l && child_l.is_leaf()) {
                attempt_add_candidate<Candidate, swap_order>(
                    query, child_l, can_collide, candidates);
            }
            if (intersects_r && child_r.is_leaf()) {
                attempt_add_candidate<Candidate, swap_order>(
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

#ifdef IPC_TOOLKIT_WITH_SIMD
    // SIMD Traversal
    // Traverses multiple queries simultaneously using SIMD.
    template <typename Candidate, bool swap_order, bool triangular>
    void traverse_lbvh_simd(
        const LBVH::Node* queries,
        const size_t first_query_leaf_idx,
        const size_t n_queries,
        const LBVH::Nodes& lbvh,
        const LBVH::RightmostLeaves& rightmost_leaves,
        const std::function<bool(size_t, size_t)>& can_collide,
        std::vector<Candidate>& candidates)
    {
        using batch_t = xs::batch<float>;
        assert(n_queries >= 1 && n_queries <= batch_t::size);

        // Load queries into single registers
        auto make_simd = [&](auto F) -> batch_t {
            // 1. Create a buffer of the correct architecture-dependent size
            alignas(xs::default_arch::alignment())
                std::array<float, batch_t::size>
                    buffer {};

#pragma unroll
            // 2. Fill the buffer, respecting the actual number of queries
            for (size_t i = 0; i < batch_t::size; ++i) {
                buffer[i] = (i < n_queries) ? F(static_cast<int>(i)) : 0.0f;
            }

            // 3. Load the buffer into the SIMD register
            return batch_t::load_aligned(buffer.data());
        };

        const auto q_min_x =
            make_simd([&](int k) { return queries[k].aabb_min.x(); });
        const auto q_min_y =
            make_simd([&](int k) { return queries[k].aabb_min.y(); });
        const auto q_min_z =
            make_simd([&](int k) { return queries[k].aabb_min.z(); });
        const auto q_max_x =
            make_simd([&](int k) { return queries[k].aabb_max.x(); });
        const auto q_max_y =
            make_simd([&](int k) { return queries[k].aabb_max.y(); });
        const auto q_max_z =
            make_simd([&](int k) { return queries[k].aabb_max.z(); });

        // Use a fixed-size array as a stack to avoid dynamic allocations
        constexpr int MAX_STACK_SIZE = 64;
        int stack[MAX_STACK_SIZE];
        int stack_ptr = 0;
        stack[stack_ptr++] = LBVH::Node::INVALID_POINTER;

        int node_idx = 0; // root
        do {
            const LBVH::Node& node = lbvh[node_idx];

            if (lbvh.size() == 1) {     // Single node case (only root)
                assert(node.is_leaf()); // Only one node, so it must be a leaf
                if constexpr (triangular) {
                    break; // No self-collision if only one node
                }
                // Check intersection with all queries simultaneously
                const xs::batch_bool<float> intersects =
                    (node.aabb_min.x() <= q_max_x)
                    & (node.aabb_min.y() <= q_max_y)
                    & (node.aabb_min.z() <= q_max_z)
                    & (q_min_x <= node.aabb_max.x())
                    & (q_min_y <= node.aabb_max.y())
                    & (q_min_z <= node.aabb_max.z());
                if (xs::any(intersects)) {
                    for (int k = 0; k < n_queries; ++k) {
                        if (intersects.get(k)) {
                            attempt_add_candidate<Candidate, swap_order>(
                                queries[k], node, can_collide, candidates);
                        }
                    }
                }
                break;
            }

            // Check left and right are valid pointers
            assert(node.is_inner());

#if defined(__GNUC__) || defined(__clang__)
            // Prefetch child nodes to reduce cache misses
            __builtin_prefetch(&lbvh[node.left], 0, 1);
            __builtin_prefetch(&lbvh[node.right], 0, 1);
#endif

            const LBVH::Node& child_l = lbvh[node.left];
            const LBVH::Node& child_r = lbvh[node.right];

            // 1. Intersect multiple queries at once
            // (child_l.min <= query.max) && (query.min <= child_l.max)
            xs::batch_bool<float> intersects_l =
                (child_l.aabb_min.x() <= q_max_x)
                & (child_l.aabb_min.y() <= q_max_y)
                & (child_l.aabb_min.z() <= q_max_z)
                & (q_min_x <= child_l.aabb_max.x())
                & (q_min_y <= child_l.aabb_max.y())
                & (q_min_z <= child_l.aabb_max.z());

            // 2. Intersect multiple queries at once
            // (child_r.min <= query.max) && (query.min <= child_r.max)
            xs::batch_bool<float> intersects_r =
                (child_r.aabb_min.x() <= q_max_x)
                & (child_r.aabb_min.y() <= q_max_y)
                & (child_r.aabb_min.z() <= q_max_z)
                & (q_min_x <= child_r.aabb_max.x())
                & (q_min_y <= child_r.aabb_max.y())
                & (q_min_z <= child_r.aabb_max.z());

            // Ignore overlap if the subtree is fully on the left-hand side
            // of all queries (triangular traversal only).
            // We use first_query_leaf_idx (the smallest query leaf index
            // in the SIMD batch) for a conservative check: if all leaves
            // in the subtree are <= the smallest query, they are also <=
            // every other query in the batch.
            if constexpr (triangular) {
                if (rightmost_leaves[node.left] <= first_query_leaf_idx) {
                    intersects_l = xs::batch_bool<float>(false);
                }
                if (rightmost_leaves[node.right] <= first_query_leaf_idx) {
                    intersects_r = xs::batch_bool<float>(false);
                }
            }

            const bool any_intersects_l = xs::any(intersects_l);
            const bool any_intersects_r = xs::any(intersects_r);

            // Query overlaps a leaf node => report collision
            if (any_intersects_l && child_l.is_leaf()) {
                for (int k = 0; k < n_queries; ++k) {
                    if constexpr (triangular) {
                        if (rightmost_leaves[node.left]
                            <= first_query_leaf_idx + k) {
                            continue;
                        }
                    }
                    if (intersects_l.get(k)) {
                        attempt_add_candidate<Candidate, swap_order>(
                            queries[k], child_l, can_collide, candidates);
                    }
                }
            }
            if (any_intersects_r && child_r.is_leaf()) {
                for (int k = 0; k < n_queries; ++k) {
                    if constexpr (triangular) {
                        if (rightmost_leaves[node.right]
                            <= first_query_leaf_idx + k) {
                            continue;
                        }
                    }
                    if (intersects_r.get(k)) {
                        attempt_add_candidate<Candidate, swap_order>(
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
        const LBVH::Nodes& source,
        const LBVH::Nodes& target,
        const LBVH::RightmostLeaves& rightmost_leaves,
        const std::function<bool(size_t, size_t)>& can_collide,
        tbb::enumerable_thread_specific<std::vector<Candidate>>& storage)
    {
#ifdef IPC_TOOLKIT_WITH_SIMD // Enable SIMD acceleration when available
        constexpr size_t SIMD_SIZE = use_simd ? xs::batch<float>::size : 1;
        static_assert(
            64 % xs::batch<float>::size == 0, "GRAIN_SIZE must be an integer");
        constexpr size_t GRAIN_SIZE =
            use_simd ? (64 / xs::batch<float>::size) : 1;
#else
        constexpr size_t SIMD_SIZE = 1;
        constexpr size_t GRAIN_SIZE = 1;
#endif

        // Calculate the offset to the first leaf node in the source BVH.
        const size_t source_leaf_offset = source.size() / 2;
        const size_t n_source_leaves = source_leaf_offset + 1;

        const size_t n_tasks =
            n_source_leaves / SIMD_SIZE + (n_source_leaves % SIMD_SIZE != 0);

        tbb::parallel_for(
            tbb::blocked_range<size_t>(size_t(0), n_tasks, GRAIN_SIZE),
            [&](const tbb::blocked_range<size_t>& r) {
                auto& local_candidates = storage.local();
#ifdef IPC_TOOLKIT_WITH_SIMD
                const size_t actual_end = // Handle tail case
                    std::min(SIMD_SIZE * r.end(), n_source_leaves);
#endif
                for (size_t i = r.begin(); i < r.end(); ++i) {
                    const size_t idx = SIMD_SIZE * i;
#ifdef IPC_TOOLKIT_WITH_SIMD
                    if constexpr (use_simd) {
                        assert(actual_end - idx >= 1);
                        traverse_lbvh_simd<Candidate, swap_order, triangular>(
                            &source[source_leaf_offset + idx], idx,
                            std::min(SIMD_SIZE, actual_end - idx), target,
                            rightmost_leaves, can_collide, local_candidates);
                    } else {
#endif
                        traverse_lbvh<Candidate, swap_order, triangular>(
                            source[source_leaf_offset + idx], idx, target,
                            rightmost_leaves, can_collide, local_candidates);
#ifdef IPC_TOOLKIT_WITH_SIMD
                    }
#endif
                }
            });
    }
} // namespace

template <typename Candidate, bool swap_order, bool triangular>
void LBVH::detect_candidates(
    const Nodes& source,
    const Nodes& target,
    const RightmostLeaves& rightmost_leaves,
    const std::function<bool(size_t, size_t)>& can_collide,
    std::vector<Candidate>& candidates)
{
    if (source.empty() || target.empty()) {
        return;
    }

    tbb::enumerable_thread_specific<std::vector<Candidate>> storage;

    independent_traversal<Candidate, swap_order, triangular>(
        source, target, rightmost_leaves, can_collide, storage);

    merge_thread_local_vectors(storage, candidates);
}

void LBVH::detect_vertex_vertex_candidates(
    std::vector<VertexVertexCandidate>& candidates) const
{
    if (vertex_bvh.size() <= 1) { // Need at least 2 vertices for a collision
        return;
    }

    IPC_TOOLKIT_PROFILE_BLOCK("LBVH::detect_vertex_vertex_candidates");

    detect_candidates(
        vertex_bvh, vertex_rightmost_leaves, can_vertices_collide, candidates);
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
        edge_bvh, vertex_bvh, vertex_rightmost_leaves,
        std::bind(&LBVH::can_edge_vertex_collide, this, _1, _2), candidates);
}

void LBVH::detect_edge_edge_candidates(
    std::vector<EdgeEdgeCandidate>& candidates) const
{
    if (edge_bvh.size() <= 1) { // Need at least 2 edges for a collision
        return;
    }

    IPC_TOOLKIT_PROFILE_BLOCK("LBVH::detect_edge_edge_candidates");

    detect_candidates(
        edge_bvh, edge_rightmost_leaves,
        std::bind(&LBVH::can_edges_collide, this, _1, _2), candidates);
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
        vertex_bvh, face_bvh, face_rightmost_leaves,
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
        face_bvh, edge_bvh, edge_rightmost_leaves,
        std::bind(&LBVH::can_edge_face_collide, this, _1, _2), candidates);
}

void LBVH::detect_face_face_candidates(
    std::vector<FaceFaceCandidate>& candidates) const
{
    if (face_bvh.size() <= 1) { // Need at least 2 faces for a collision
        return;
    }

    IPC_TOOLKIT_PROFILE_BLOCK("LBVH::detect_face_face_candidates");
    detect_candidates(
        face_bvh, face_rightmost_leaves,
        std::bind(&LBVH::can_faces_collide, this, _1, _2), candidates);
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