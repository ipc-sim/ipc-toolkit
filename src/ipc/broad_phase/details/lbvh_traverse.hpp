#pragma once

#include <ipc/config.hpp>
#include <ipc/broad_phase/lbvh.hpp>

#include <cassert>
#include <cstdint>

namespace ipc::details {

/// @brief Descends a target BVH for one query leaf, reporting every target leaf
/// whose AABB overlaps the query.
///
/// A stackless-style descent with an explicit stack: at each inner node the
/// overlapping children are handled immediately if they are leaves, descended
/// into if only one is inner, and the right one postponed on the stack if both
/// are. The root lives at index 0, and LBVH::Node::INVALID_POINTER (which is
/// also 0) doubles as the stack's bottom sentinel -- popping it ends the walk,
/// because no node other than the root ever lives at index 0.
///
/// This is shared by ipc::LBVH and ipc::cuda::LBVH. What differs between them
/// is only what happens on an overlap, which is why that is a policy: the host
/// filters and appends to a std::vector, while the device filters against the
/// mesh connectivity and appends through an atomic counter.
///
/// @tparam triangular Self-collision: skip any subtree lying entirely to the
/// left of the query, so each unordered pair is reported exactly once.
/// @tparam Emit Callable (const LBVH::Node& leaf) -> void, invoked for every
/// overlapping target leaf. It owns both the collision filtering and the
/// recording of the pair.
///
/// @param query The querying leaf node.
/// @param query_leaf_idx The query's position in its own Morton-sorted leaf
/// order. Used only by the triangular skip.
/// @param target The target BVH's nodes, root at index 0.
/// @param target_size The number of nodes in the target BVH.
/// @param target_rightmost The target's per-node rightmost-leaf indices. Used
/// only by the triangular skip.
/// @param emit The per-overlap callback.
template <bool triangular, typename Emit>
IPC_TOOLKIT_HOST_DEVICE void traverse_lbvh(
    const LBVH::Node& query,
    const int query_leaf_idx,
    const LBVH::Node* target,
    const int target_size,
    const int32_t* target_rightmost,
    Emit&& emit)
{
    // A fixed-size stack keeps the descent free of dynamic allocation.
    constexpr int MAX_STACK_SIZE = 64;
    int stack[MAX_STACK_SIZE];
    int stack_ptr = 0;
    stack[stack_ptr++] = LBVH::Node::INVALID_POINTER;

    int node_idx = 0; // root
    do {
        const LBVH::Node& node = target[node_idx];

        if (target_size == 1) { // only the root, which is therefore a leaf
            assert(node.is_leaf());
            if constexpr (triangular) {
                break; // a lone primitive cannot collide with itself
            }
            if (node.intersects(query)) {
                emit(node);
            }
            break;
        }

        assert(node.is_inner()); // so .left and .right are valid pointers

#if !defined(__CUDA_ARCH__) && (defined(__GNUC__) || defined(__clang__))
        // Prefetch the children to reduce cache misses. The device needs no
        // equivalent; it hides the latency with its other resident warps.
        __builtin_prefetch(&target[node.left], 0, 1);
        __builtin_prefetch(&target[node.right], 0, 1);
#endif

        const LBVH::Node& child_l = target[node.left];
        const LBVH::Node& child_r = target[node.right];
        bool intersects_l = child_l.intersects(query);
        bool intersects_r = child_r.intersects(query);

        // Ignore a subtree lying entirely to the query's left; that pair is
        // reported when the other primitive is the query instead.
        if constexpr (triangular) {
            if (intersects_l && target_rightmost[node.left] <= query_leaf_idx) {
                intersects_l = false;
            }
            if (intersects_r
                && target_rightmost[node.right] <= query_leaf_idx) {
                intersects_r = false;
            }
        }

        const bool l_leaf = child_l.is_leaf();
        const bool r_leaf = child_r.is_leaf();

        // An overlapped leaf is a candidate.
        if (intersects_l && l_leaf) {
            emit(child_l);
        }
        if (intersects_r && r_leaf) {
            emit(child_r);
        }

        // An overlapped inner node is descended into.
        const bool traverse_l = intersects_l && !l_leaf;
        const bool traverse_r = intersects_r && !r_leaf;

        if (!traverse_l && !traverse_r) {
            assert(stack_ptr > 0);
            node_idx = stack[--stack_ptr];
        } else {
            node_idx = traverse_l ? node.left : node.right;
            if (traverse_l && traverse_r) {
                assert(stack_ptr < MAX_STACK_SIZE);
                stack[stack_ptr++] = node.right; // postpone the right child
            }
        }
    } while (node_idx != LBVH::Node::INVALID_POINTER);
}

} // namespace ipc::details
