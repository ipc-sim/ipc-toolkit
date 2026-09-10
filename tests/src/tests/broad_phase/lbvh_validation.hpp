#pragma once

// Structural validators for ipc::LBVH node arrays, shared by the CPU tests
// (test_lbvh.cpp) and the GPU parity tests (test_gpu_lbvh.cu) so both trees
// are held to the same predicate. Host-only: no CUDA here.

#include <ipc/broad_phase/lbvh.hpp>

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cstdint>
#include <vector>

namespace ipc::tests {

/// @brief Whether a parent's AABB is exactly the union of its children's.
///
/// Exact, not approximate: the build computes the parent as the component-wise
/// float min/max of its children on the host and the device alike, so anything
/// short of bit equality means a parent was combined from the wrong children or
/// was written before both had arrived.
inline bool is_aabb_union(
    const LBVH::Node& parent,
    const LBVH::Node& child_a,
    const LBVH::Node& child_b)
{
    return (parent.aabb_min == child_a.aabb_min.min(child_b.aabb_min)).all()
        && (parent.aabb_max == child_a.aabb_max.max(child_b.aabb_max)).all();
}

/// @brief Recursively verify that every node below @p index is reached exactly
/// once and that every internal node's AABB is the union of its children's,
/// collecting the primitive ids of the leaves reached.
///
/// A revisit is fatal (REQUIRE), not just a failure: a shared child would
/// otherwise be re-traversed exponentially, and a cycle forever.
inline void traverse_lbvh_nodes(
    const LBVH::Nodes& nodes,
    const int32_t index,
    std::vector<bool>& visited,
    std::vector<int32_t>& reached_leaves)
{
    REQUIRE(index >= 0);
    REQUIRE(index < int32_t(nodes.size()));
    const LBVH::Node& node = nodes[index];
    CHECK(node.is_valid());
    REQUIRE(!visited[index]);
    visited[index] = true;

    if (node.is_leaf()) {
        reached_leaves.push_back(node.primitive_id);
        return;
    }

    const LBVH::Node& child_a = nodes[node.left];
    const LBVH::Node& child_b = nodes[node.right];
    {
        CAPTURE(
            index, node.left, node.right, node.aabb_min.transpose(),
            child_a.aabb_min.transpose(), child_b.aabb_min.transpose(),
            node.aabb_max.transpose(), child_a.aabb_max.transpose(),
            child_b.aabb_max.transpose());
        CHECK(is_aabb_union(node, child_a, child_b));
    }
    traverse_lbvh_nodes(nodes, node.left, visited, reached_leaves);
    traverse_lbvh_nodes(nodes, node.right, visited, reached_leaves);
}

/// @brief Validate a built LBVH: 2n - 1 nodes for n leaves, every node
/// reachable exactly once from the root at index 0, every internal AABB the
/// union of its children's, and the leaves holding exactly the primitive ids
/// {0, ..., n - 1}. A single-node tree (one primitive) is valid.
inline void check_valid_lbvh_nodes(const LBVH::Nodes& nodes)
{
    REQUIRE(!nodes.empty());
    REQUIRE(nodes.size() % 2 == 1);
    const size_t n_leaves = (nodes.size() + 1) / 2;

    std::vector<bool> visited(nodes.size(), false);
    std::vector<int32_t> reached_leaves;
    traverse_lbvh_nodes(nodes, 0, visited, reached_leaves);
    CHECK(
        std::all_of(visited.begin(), visited.end(), [](bool v) { return v; }));

    REQUIRE(reached_leaves.size() == n_leaves);
    std::sort(reached_leaves.begin(), reached_leaves.end());
    for (size_t i = 0; i < reached_leaves.size(); ++i) {
        CHECK(reached_leaves[i] == int32_t(i));
    }
}

/// @brief Validate a tree built by one implementation against the same tree
/// built by another from the same boxes: it must be a valid LBVH of the same
/// size, and its root AABB -- an order-independent union of identically
/// inflated boxes -- must be exactly equal.
///
/// The size check comes first, so an implementation that returns no nodes
/// fails rather than passing vacuously.
inline void
check_lbvh_nodes_match(const LBVH::Nodes& nodes, const LBVH::Nodes& reference)
{
    REQUIRE(nodes.size() == reference.size());
    check_valid_lbvh_nodes(nodes);
    CHECK((nodes[0].aabb_min == reference[0].aabb_min).all());
    CHECK((nodes[0].aabb_max == reference[0].aabb_max).all());
}

} // namespace ipc::tests
