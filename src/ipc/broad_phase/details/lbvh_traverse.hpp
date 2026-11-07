#pragma once

#include <ipc/config.hpp>
#include <ipc/broad_phase/lbvh.hpp>
#include <ipc/math/morton.hpp>  // for MORTON_KEY_BITS
#include <ipc/utils/logger.hpp> // for log_and_throw_error (host only)
#include <ipc/utils/simd.hpp>   // for any(bool)

#include <cassert>
#include <cstdint>
#include <type_traits>

namespace ipc::details {

/// @brief Descends a target BVH for one query -- or one batch of queries --
/// reporting every target leaf whose AABB overlaps.
///
/// A stackless-style descent with an explicit stack: at each inner node the
/// overlapping children are handled immediately if they are leaves, descended
/// into if only one is inner, and the right one postponed on the stack if both
/// are. The root lives at index 0, and LBVH::Node::INVALID_POINTER (which is
/// also 0) doubles as the stack's bottom sentinel -- popping it ends the walk,
/// because no node other than the root ever lives at index 0.
///
/// This is shared by ipc::LBVH -- both its scalar and its SIMD traversal -- and
/// ipc::cuda::LBVH. What differs between them is how a node is tested against
/// the query and what happens on an overlap, which is why both are policies:
/// the scalar host and the device test one query and get a bool; the SIMD host
/// tests a batch of queries at once and gets a lane mask. On an overlap the
/// host filters and appends to a std::vector, while the device filters against
/// the mesh connectivity and appends through an atomic counter.
///
/// @tparam triangular Self-collision: skip any subtree lying entirely to the
/// left of the query, so each unordered pair is reported exactly once.
/// @tparam Intersects Callable (const LBVH::Node&) -> Mask, where Mask is bool
/// for a single query or a lane mask for a batch of queries. Either must
/// support any(mask) -- ipc::any for bool, xsimd::any by ADL for a batch --
/// and Mask(false).
/// @tparam Emit Callable (const LBVH::Node& leaf, int leaf_idx, const Mask&)
/// -> void, invoked for every target leaf overlapping at least one query. It
/// owns both the collision filtering and the recording of the pair. For a
/// batch, it must re-check the triangular skip per lane against
/// target_rightmost[leaf_idx], because the skip below is conservative for the
/// batch as a whole.
///
/// @param query_leaf_idx The query's position in its own Morton-sorted leaf
/// order -- for a batch, the smallest position in it. Used only by the
/// triangular skip.
/// @param target The target BVH's nodes, root at index 0.
/// @param target_size The number of nodes in the target BVH.
/// @param target_rightmost The target's per-node rightmost-leaf indices. Used
/// only by the triangular skip.
/// @param intersects The per-node overlap test.
/// @param emit The per-overlap callback.
template <bool triangular, typename Intersects, typename Emit>
IPC_TOOLKIT_HOST_DEVICE void traverse_lbvh(
    const int query_leaf_idx,
    const LBVH::Node* target,
    const int target_size,
    const int32_t* target_rightmost,
    Intersects&& intersects,
    Emit&& emit)
{
    using Mask = std::decay_t<decltype(intersects(*target))>;

    // A fixed-size stack keeps the descent free of dynamic allocation. Every
    // internal node on a root-to- leaf path splits at a distinct, strictly
    // increasing prefix length of the MORTON_KEY_BITS-bit key, so a path holds
    // at most MORTON_KEY_BITS internal nodes, and at most one right child per
    // internal node on the current path is ever pending -- plus the sentinel.
    // The overflow check below can therefore never fire on a well-formed tree;
    // it turns a malformed one into a hard failure instead of an out-of-bounds
    // write.
    constexpr int MAX_STACK_SIZE = MORTON_KEY_BITS + 1;
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
            const Mask mask = intersects(node);
            if (any(mask)) {
                emit(node, node_idx, mask);
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
        Mask intersects_l = intersects(child_l);
        Mask intersects_r = intersects(child_r);

        // Ignore a subtree lying entirely to the query's left; that pair is
        // reported when the other primitive is the query instead. For a batch
        // this uses its smallest query position, so it is conservative: a
        // subtree left of every query is skipped here, and emit() re-checks
        // the rest per lane.
        if constexpr (triangular) {
            if constexpr (std::is_same_v<Mask, bool>) {
                // A scalar overlap result is a free branch, so test it first
                // and skip the rightmost-leaf load when nothing overlaps.
                if (intersects_l
                    && target_rightmost[node.left] <= query_leaf_idx) {
                    intersects_l = false;
                }
                if (intersects_r
                    && target_rightmost[node.right] <= query_leaf_idx) {
                    intersects_r = false;
                }
            } else {
                // For a batch, any(mask) is a movemask plus a second
                // data-dependent branch per child, which costs more than the
                // one load it would save -- so test the position alone.
                if (target_rightmost[node.left] <= query_leaf_idx) {
                    intersects_l = Mask(false);
                }
                if (target_rightmost[node.right] <= query_leaf_idx) {
                    intersects_r = Mask(false);
                }
            }
        }

        const bool any_l = any(intersects_l);
        const bool any_r = any(intersects_r);

        // An overlapped leaf is a candidate.
        if (any_l && child_l.is_leaf()) {
            emit(child_l, node.left, intersects_l);
        }
        if (any_r && child_r.is_leaf()) {
            emit(child_r, node.right, intersects_r);
        }

        // An overlapped inner node is descended into.
        const bool traverse_l = any_l && !child_l.is_leaf();
        const bool traverse_r = any_r && !child_r.is_leaf();

        if (!traverse_l && !traverse_r) {
            assert(stack_ptr > 0);
            node_idx = stack[--stack_ptr];
        } else {
            node_idx = traverse_l ? node.left : node.right;
            if (traverse_l && traverse_r) {
                if (stack_ptr >= MAX_STACK_SIZE) {
                    // Unreachable on a well-formed tree (see MAX_STACK_SIZE);
                    // fail hard rather than write past the stack.
#ifdef __CUDA_ARCH__
                    __trap();
#else
                    log_and_throw_error(
                        "ipc::details::traverse_lbvh: traversal stack "
                        "overflow; the BVH is deeper than the Morton key "
                        "width allows and so is malformed");
#endif
                }
                stack[stack_ptr++] = node.right; // postpone the right child
            }
        }
    } while (node_idx != LBVH::Node::INVALID_POINTER);
}

} // namespace ipc::details
