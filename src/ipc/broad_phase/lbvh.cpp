#include "lbvh.hpp"

#include <ipc/broad_phase/details/connectivity_filters.hpp>
#include <ipc/broad_phase/details/lbvh_build.hpp>
#include <ipc/broad_phase/details/lbvh_traverse.hpp>
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

        const Eigen::Array3d mesh_width_inv =
            morton_domain_width_inv(mesh_aabb.min, mesh_aabb.max);
        tbb::parallel_for(size_t(0), boxes.size(), [&](size_t i) {
            const auto& box = boxes[i];

            morton_codes[i].morton_code = morton_code(
                0.5 * (box.min + box.max), mesh_aabb.min, mesh_width_inv, dim);
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
    // hierarchy and computes bounding boxes. See
    // ipc::details::build_hierarchy_from_leaf(), shared with the device build
    // in ipc::cuda::LBVH.
    std::atomic<int> root_idx(-1);
    {
        IPC_TOOLKIT_PROFILE_BLOCK("build_hierarchy_and_boxes");
        tbb::parallel_for(0, N_LEAVES, [&](int i) {
            const size_t box_id = morton_codes[i].box_id;
            details::init_leaf_node(
                i, N_LEAVES, box_id, boxes[box_id].min, boxes[box_id].max,
                lbvh.data(), rightmost_leaves.data());

            const int root = details::build_hierarchy_from_leaf(
                i, N_LEAVES, [&](int k) { return morton_codes[k].morton_code; },
                lbvh.data(), rightmost_leaves.data(), construction_infos.data(),
                // std::atomic's post-increment is sequentially consistent, so
                // it already orders this thread's writes before the arrival
                // and the arrival before its later reads.
                [](std::atomic<int>& count) { return count++; });

            if (root >= 0) {
                // Only one thread should ever reach the root.
                int expected = -1;
                [[maybe_unused]] const bool set =
                    root_idx.compare_exchange_strong(expected, root);
                assert(set);
            }
        });
    }

    // --- Move the root to index 0 so traversal can start there. ---
    // In the Apetrei layout the root's index equals the global split position,
    // which is generally != 0.
    const int root = root_idx.load();
    if (root > 0) {
        IPC_TOOLKIT_PROFILE_BLOCK("swap_root_to_zero");
        details::swap_root_to_zero(lbvh.data(), rightmost_leaves.data(), root);

        tbb::parallel_for(size_t(0), lbvh.size(), [&](size_t i) {
            details::patch_left_pointer(lbvh[i], root);
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

    /// Scalar traversal: descend the target BVH for one query leaf and record
    /// every overlapping, filter-passing pair. The descent itself is
    /// ipc::details::traverse_lbvh(), shared with the SIMD traversal below and
    /// with ipc::cuda::LBVH; only the overlap test and the emission are ours.
    template <typename Candidate, bool swap_order, bool triangular>
    void traverse_lbvh(
        const LBVH::Node& query,
        const size_t query_leaf_idx,
        const LBVH::Nodes& lbvh,
        const LBVH::RightmostLeaves& rightmost_leaves,
        const std::function<bool(size_t, size_t)>& can_collide,
        std::vector<Candidate>& candidates)
    {
        details::traverse_lbvh<triangular>(
            int(query_leaf_idx), lbvh.data(), int(lbvh.size()),
            rightmost_leaves.data(),
            [&](const LBVH::Node& node) { return node.intersects(query); },
            [&](const LBVH::Node& leaf, const int /*leaf_idx*/,
                const bool /*intersects*/) {
                attempt_add_candidate<Candidate, swap_order>(
                    query, leaf, can_collide, candidates);
            });
    }

#ifdef IPC_TOOLKIT_WITH_SIMD
    /// SIMD traversal: descend the target BVH once for a batch of up to
    /// xs::batch<float>::size query leaves, testing every node against all of
    /// them at once. The descent is the same ipc::details::traverse_lbvh() as
    /// the scalar path; the overlap test returns a lane mask instead of a bool,
    /// and the emission fans a leaf out to the lanes that overlap it.
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
        using mask_t = xs::batch_bool<float>;
        assert(n_queries >= 1 && n_queries <= batch_t::size);

        // Load queries into single registers
        auto make_simd = [&](auto F) -> batch_t {
            // 1. Create a buffer of the correct architecture-dependent size
            alignas(xs::default_arch::alignment())
                std::array<float, batch_t::size>
                    buffer {};

#if defined(__clang__)
#pragma unroll
#elif defined(__GNUC__)
#pragma GCC unroll 16
#endif
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

        details::traverse_lbvh<triangular>(
            int(first_query_leaf_idx), lbvh.data(), int(lbvh.size()),
            rightmost_leaves.data(),
            // Intersect all queries at once:
            // (node.min <= query.max) && (query.min <= node.max)
            [&](const LBVH::Node& node) -> mask_t {
                return (node.aabb_min.x() <= q_max_x)
                    & (node.aabb_min.y() <= q_max_y)
                    & (node.aabb_min.z() <= q_max_z)
                    & (q_min_x <= node.aabb_max.x())
                    & (q_min_y <= node.aabb_max.y())
                    & (q_min_z <= node.aabb_max.z());
            },
            [&](const LBVH::Node& leaf, const int leaf_idx,
                const mask_t& intersects) {
                for (size_t k = 0; k < n_queries; ++k) {
                    if constexpr (triangular) {
                        // The shared descent skipped subtrees left of the
                        // batch's FIRST query; finish the check per lane.
                        if (rightmost_leaves[leaf_idx]
                            <= first_query_leaf_idx + k) {
                            continue;
                        }
                    }
                    if (intersects.get(k)) {
                        attempt_add_candidate<Candidate, swap_order>(
                            queries[k], leaf, can_collide, candidates);
                    }
                }
            });
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
    candidates.clear();
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
    candidates.clear();
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
    candidates.clear();
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
    candidates.clear();
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
    candidates.clear();
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
    candidates.clear();
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

    return details::can_edge_vertex_collide(e0i, e1i, vi, can_vertices_collide);
}

bool LBVH::can_edges_collide(size_t eai, size_t ebi) const
{
    assert(eai < edge_vertex_ids.size());
    const auto& [ea0i, ea1i] = edge_vertex_ids[eai];
    assert(ebi < edge_vertex_ids.size());
    const auto& [eb0i, eb1i] = edge_vertex_ids[ebi];

    return details::can_edges_collide(
        ea0i, ea1i, eb0i, eb1i, can_vertices_collide);
}

bool LBVH::can_face_vertex_collide(size_t fi, size_t vi) const
{
    assert(fi < face_vertex_ids.size());
    const auto& [f0i, f1i, f2i] = face_vertex_ids[fi];

    return details::can_face_vertex_collide(
        f0i, f1i, f2i, vi, can_vertices_collide);
}

bool LBVH::can_edge_face_collide(size_t ei, size_t fi) const
{
    assert(ei < edge_vertex_ids.size());
    const auto& [e0i, e1i] = edge_vertex_ids[ei];
    assert(fi < face_vertex_ids.size());
    const auto& [f0i, f1i, f2i] = face_vertex_ids[fi];

    return details::can_edge_face_collide(
        e0i, e1i, f0i, f1i, f2i, can_vertices_collide);
}

bool LBVH::can_faces_collide(size_t fai, size_t fbi) const
{
    assert(fai < face_vertex_ids.size());
    const auto& [fa0i, fa1i, fa2i] = face_vertex_ids[fai];
    assert(fbi < face_vertex_ids.size());
    const auto& [fb0i, fb1i, fb2i] = face_vertex_ids[fbi];

    return details::can_faces_collide(
        fa0i, fa1i, fa2i, fb0i, fb1i, fb2i, can_vertices_collide);
}

} // namespace ipc