#pragma once

#include <ipc/config.hpp>
#include <ipc/broad_phase/lbvh.hpp>
#include <ipc/math/morton.hpp>

#include <Eigen/Core>

#include <cassert>
#include <cmath>
#include <cstdint>

namespace ipc::details {

// The LBVH construction of Apetrei [2014], shared by ipc::LBVH and
// ipc::cuda::LBVH.
//
// Everything here is host/device and addresses the tree through raw pointers,
// so the CPU can drive it from a tbb::parallel_for over std::vectors and the
// GPU from a kernel over thrust::device_vectors. What each platform still owns
// is only the parallel launch, the sort, and the two policies these take: how
// to read a sorted Morton code, and how to perform the atomic arrival.

/// @brief Rounds a double AABB outward to the smallest enclosing float AABB.
///
/// Each corner is nudged to the next representable float away from the box, so
/// the float AABB always encloses the double one and never clips a primitive.
///
/// @param box_min The minimum corner of the double AABB.
/// @param box_max The maximum corner of the double AABB.
/// @param[out] node The node whose AABB is set.
IPC_TOOLKIT_HOST_DEVICE inline void set_inflated_aabb(
    const Eigen::Array3d& box_min,
    const Eigen::Array3d& box_max,
    LBVH::Node& node)
{
    for (int k = 0; k < 3; ++k) {
        node.aabb_min[k] =
            nextafterf(static_cast<float>(box_min[k]), -INFINITY);
        node.aabb_max[k] = nextafterf(static_cast<float>(box_max[k]), INFINITY);
    }
}

/// @brief Initializes the leaf node for sorted position i.
///
/// Leaves occupy the upper half of the node array, at n_leaves - 1 + i.
///
/// @param i The leaf's position in the Morton-sorted order.
/// @param n_leaves The number of leaves.
/// @param box_id The id of the primitive this leaf holds.
/// @param box_min The minimum corner of the primitive's AABB.
/// @param box_max The maximum corner of the primitive's AABB.
/// @param[out] nodes The BVH nodes.
/// @param[out] rightmost_leaves The per-node rightmost-leaf indices.
IPC_TOOLKIT_HOST_DEVICE inline void init_leaf_node(
    const int i,
    const int n_leaves,
    const index_t box_id,
    const Eigen::Array3d& box_min,
    const Eigen::Array3d& box_max,
    LBVH::Node* nodes,
    int32_t* rightmost_leaves)
{
    LBVH::Node leaf;
    set_inflated_aabb(box_min, box_max, leaf);
    leaf.primitive_id = static_cast<int32_t>(box_id);
    leaf.is_inner_marker = 0;

    const int leaf_idx = n_leaves - 1 + i;
    nodes[leaf_idx] = leaf;
    rightmost_leaves[leaf_idx] = i; // a leaf's rightmost leaf is itself
}

/// @brief Returns the length of the common leading-bit prefix of the sorted
/// Morton codes at positions i and j, or -1 when j is out of bounds.
/// @tparam CodeAt Callable (int) -> uint64_t returning a sorted Morton code.
/// @param code_at The sorted Morton code accessor.
/// @param n_leaves The number of codes.
/// @param i The first sorted position.
/// @param j The second sorted position.
/// @return The common-prefix length, or -1 when j is out of bounds.
template <typename CodeAt>
IPC_TOOLKIT_HOST_DEVICE inline int
delta(CodeAt&& code_at, const int n_leaves, const int i, const int j)
{
    if (j < 0 || j >= n_leaves) {
        return -1;
    }
    return morton_common_prefix(code_at(i), i, code_at(j), j);
}

/// @brief Whether the subtree spanning [left_key, right_key] is its parent's
/// left child.
///
/// The two candidate parents are internal node right_key (which would make
/// this the left child) and internal node left_key - 1 (the right child).
/// delta() grows with the similarity of the codes, so the nearer ancestor is
/// the one with the LARGER delta -- hence ">".
///
/// At the boundaries only one candidate exists: a range starting at 0 has no
/// node -1 to its left, and a range ending at n_leaves - 1 has no node
/// n_leaves - 1 to its right.
///
/// @tparam CodeAt Callable (int) -> uint64_t returning a sorted Morton code.
/// @param code_at The sorted Morton code accessor.
/// @param n_leaves The number of leaves.
/// @param left_key The left endpoint of the subtree's sorted-key range.
/// @param right_key The right endpoint of the subtree's sorted-key range.
/// @return Whether this subtree is its parent's left child.
template <typename CodeAt>
IPC_TOOLKIT_HOST_DEVICE inline bool is_left_child(
    CodeAt&& code_at,
    const int n_leaves,
    const int left_key,
    const int right_key)
{
    return left_key == 0
        || (right_key != n_leaves - 1
            && delta(code_at, n_leaves, right_key, right_key + 1)
                > delta(code_at, n_leaves, left_key - 1, left_key));
}

/// @brief Walks one leaf up to the root, building the hierarchy, the internal
/// AABBs and the rightmost-leaf indices (Apetrei [2014], Fig. 2).
///
/// Each leaf's thread climbs toward the root, choosing its parent in O(1) from
/// the delta values at the two ends of its current key range. At every parent
/// the first of the two arriving threads stops and the second continues, so
/// whoever continues knows both children are complete.
///
/// In this layout internal node j always splits between sorted keys j and
/// j + 1, so the root is generally NOT at index 0; swap_root_to_zero() moves
/// it there afterwards, which is what the traversal expects.
///
/// @tparam Counter The visitation counter's type (see
/// ipc::LBVH::ConstructionInfo).
/// @tparam CodeAt Callable (int) -> uint64_t returning a sorted Morton code.
/// @tparam Arrive Callable (Counter&) -> int that atomically increments the
/// counter and returns its previous value. It must order this thread's earlier
/// writes before the increment, and -- when it returns nonzero, so this thread
/// continues -- order the increment before this thread's later reads. Without
/// both halves the continuing thread can read a stale sibling.
///
/// @param i The leaf's position in the Morton-sorted order.
/// @param n_leaves The number of leaves.
/// @param code_at The sorted Morton code accessor.
/// @param[in,out] nodes The BVH nodes; the leaves must already be initialized.
/// @param[in,out] rightmost_leaves The per-node rightmost-leaf indices.
/// @param[in,out] infos The per-node construction scratch, zero-initialized.
/// @param arrive The atomic arrival gate.
/// @return The root's index if this leaf's walk reached the root, else -1.
template <typename Counter, typename CodeAt, typename Arrive>
IPC_TOOLKIT_HOST_DEVICE int build_hierarchy_from_leaf(
    const int i,
    const int n_leaves,
    CodeAt&& code_at,
    LBVH::Node* nodes,
    int32_t* rightmost_leaves,
    LBVH::ConstructionInfo<Counter>* infos,
    Arrive&& arrive)
{
    // A single-leaf tree is its own root and has no internal nodes to build.
    if (n_leaves == 1) {
        return i == 0 ? 0 : -1;
    }

    // Invariant: the current subtree covers the sorted-key range
    // [left_key, right_key].
    int left_key = i;
    int right_key = i;
    int current_node = n_leaves - 1 + i;

    while (true) {
        const bool is_child_a =
            is_left_child(code_at, n_leaves, left_key, right_key);
        const int parent = is_child_a ? right_key : left_key - 1;

        // Write the child pointer and the range endpoint onto the parent. The
        // left child writes .left and the left endpoint, the right child
        // writes .right and the right endpoint.
        if (is_child_a) {
            nodes[parent].left = current_node;
            infos[parent].left_range = left_key;
        } else {
            nodes[parent].right = current_node;
            infos[parent].right_range = right_key;
        }

        if (arrive(infos[parent].visitation_count) == 0) {
            return -1; // first thread to arrive here -> done
        }

        // Second thread to arrive: both children are complete, so their AABBs
        // and rightmost leaves can be combined into the parent's.
        assert(nodes[parent].is_inner());
        const LBVH::Node& child_a = nodes[nodes[parent].left];
        const LBVH::Node& child_b = nodes[nodes[parent].right];
        nodes[parent].aabb_min = child_a.aabb_min.min(child_b.aabb_min);
        nodes[parent].aabb_max = child_a.aabb_max.max(child_b.aabb_max);

        const int32_t rightmost_a = rightmost_leaves[nodes[parent].left];
        const int32_t rightmost_b = rightmost_leaves[nodes[parent].right];
        rightmost_leaves[parent] =
            rightmost_a > rightmost_b ? rightmost_a : rightmost_b;

        // Reconstruct the parent's full key range and continue upward.
        left_key = infos[parent].left_range;
        right_key = infos[parent].right_range;
        current_node = parent;

        if (left_key == 0 && right_key == n_leaves - 1) {
            return current_node; // the root's AABB is complete
        }
    }
}

/// @brief Swaps the node and rightmost-leaf entries at index 0 and the root, so
/// that traversal can start at index 0.
///
/// The root is never any node's child, so no pointer needs rewriting to reach
/// its new home at 0. Pointers that referenced the old node 0 do, which is what
/// patch_left_pointer() handles.
///
/// @param[in,out] nodes The BVH nodes.
/// @param[in,out] rightmost_leaves The per-node rightmost-leaf indices.
/// @param root The root's index, which must be greater than 0.
IPC_TOOLKIT_HOST_DEVICE inline void
swap_root_to_zero(LBVH::Node* nodes, int32_t* rightmost_leaves, const int root)
{
    assert(root > 0);

    const LBVH::Node node = nodes[0];
    nodes[0] = nodes[root];
    nodes[root] = node;

    const int32_t rightmost = rightmost_leaves[0];
    rightmost_leaves[0] = rightmost_leaves[root];
    rightmost_leaves[root] = rightmost;
}

/// @brief Rewrites a left pointer that referenced the old node 0 to the root's
/// new location, after swap_root_to_zero().
///
/// Apetrei's layout guarantees node 0's subtree has left_key == 0, so node 0 is
/// only ever written as a LEFT child. Two things follow: only .left pointers
/// need patching, and swapping node 0 away cannot leave a node whose .right --
/// which aliases is_inner_marker -- is 0 and so reads as a leaf.
///
/// @param[in,out] node The node to patch.
/// @param root The new location of the old node 0.
IPC_TOOLKIT_HOST_DEVICE inline void
patch_left_pointer(LBVH::Node& node, const int root)
{
    if (node.is_inner() && node.left == 0) {
        node.left = root;
    }
}

} // namespace ipc::details
