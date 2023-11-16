#include "hash_grid.hpp"

#include <ipc/broad_phase/voxel_size_heuristic.hpp>
#include <ipc/utils/merge_thread_local.hpp>
#include <ipc/utils/logger.hpp>

#include <tbb/enumerable_thread_specific.h>
#include <tbb/blocked_range2d.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>

#include <algorithm> // std::min/max

#define IPC_TOOLKIT_HASH_GRID_USE_SORT_UNIQUE // else use unordered_set

namespace ipc {

void HashGrid::build(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const double inflation_radius)
{
    BroadPhase::build(vertices, edges, faces, inflation_radius);
    // BroadPhase::build also calls clear()

    ArrayMax3d mesh_min = vertices.colwise().minCoeff().array();
    ArrayMax3d mesh_max = vertices.colwise().maxCoeff().array();
    AABB::conservative_inflation(mesh_min, mesh_max, inflation_radius);

    double cell_size =
        suggest_good_voxel_size(vertices, edges, inflation_radius);
    assert(std::isfinite(cell_size));
    resize(mesh_min, mesh_max, cell_size);

    insert_boxes();
}

void HashGrid::build(
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const double inflation_radius)
{
    BroadPhase::build(vertices_t0, vertices_t1, edges, faces, inflation_radius);
    // BroadPhase::build also calls clear()

    const ArrayMax3d mesh_min_t0 = vertices_t0.colwise().minCoeff();
    const ArrayMax3d mesh_max_t0 = vertices_t0.colwise().maxCoeff();
    const ArrayMax3d mesh_min_t1 = vertices_t1.colwise().minCoeff();
    const ArrayMax3d mesh_max_t1 = vertices_t1.colwise().maxCoeff();

    ArrayMax3d mesh_min = mesh_min_t0.min(mesh_min_t1);
    ArrayMax3d mesh_max = mesh_max_t0.max(mesh_max_t1);
    AABB::conservative_inflation(mesh_min, mesh_max, inflation_radius);

    double cell_size = suggest_good_voxel_size(
        vertices_t0, vertices_t1, edges, inflation_radius);
    assert(std::isfinite(cell_size));
    resize(mesh_min, mesh_max, cell_size);

    insert_boxes();
}

void HashGrid::resize(
    const ArrayMax3d& min, const ArrayMax3d& max, double cellSize)
{
    assert(cellSize != 0.0);
    assert(std::isfinite(cellSize));
    m_cellSize = cellSize;
    m_domainMin = min;
    m_domainMax = max;
    m_gridSize = ((max - min) / m_cellSize).ceil().cast<int>().max(1);
    logger().trace(
        "hash-grid resized with a size of {:d}x{:d}x{:d}", m_gridSize[0],
        m_gridSize[1], m_gridSize.size() == 3 ? m_gridSize[2] : 1);
}

void HashGrid::insert_boxes()
{
    insert_boxes(this->vertex_boxes, vertex_items);
    insert_boxes(this->edge_boxes, edge_items);
    insert_boxes(this->face_boxes, face_items);
}

void HashGrid::insert_boxes(
    const std::vector<AABB>& boxes, std::vector<HashItem>& items) const
{
    tbb::enumerable_thread_specific<std::vector<HashItem>> storage;

    tbb::parallel_for(
        tbb::blocked_range<long>(0l, long(boxes.size())),
        [&](const tbb::blocked_range<long>& range) {
            auto& local_items = storage.local();
            for (long i = range.begin(); i != range.end(); i++) {
                insert_box(boxes[i], i, local_items);
            }
        });

    merge_thread_local_vectors(storage, items);

    // Sorted all they (key, value) pairs, where key is the hash key, and
    // value is the element index
    tbb::parallel_sort(items.begin(), items.end());
}

void HashGrid::insert_box(
    const AABB& aabb, const long id, std::vector<HashItem>& items) const
{
    ArrayMax3i int_min = ((aabb.min - m_domainMin) / m_cellSize).cast<int>();
    // We can round down to -1, but not less
    assert((int_min >= -1).all());
    assert((int_min <= m_gridSize).all());
    int_min = int_min.max(0).min(m_gridSize - 1);

    ArrayMax3i int_max = ((aabb.max - m_domainMin) / m_cellSize).cast<int>();
    assert((int_max >= -1).all());
    assert((int_max <= m_gridSize).all());
    int_max = int_max.max(0).min(m_gridSize - 1);
    assert((int_min <= int_max).all());

    int min_z = int_min.size() == 3 ? int_min.z() : 0;
    int max_z = int_max.size() == 3 ? int_max.z() : 0;
    for (int x = int_min.x(); x <= int_max.x(); ++x) {
        for (int y = int_min.y(); y <= int_max.y(); ++y) {
            for (int z = min_z; z <= max_z; ++z) {
                items.emplace_back(hash(x, y, z), id);
            }
        }
    }
}

template <typename Candidate>
void HashGrid::detect_candidates(
    const std::vector<HashItem>& items0,
    const std::vector<HashItem>& items1,
    const std::vector<AABB>& boxes0,
    const std::vector<AABB>& boxes1,
    const std::function<bool(size_t, size_t)>& can_collide,
    std::vector<Candidate>& candidates) const
{
    // Entries with the same key means they share a cell (that cell index
    // hashes to the same key) and should be flagged for low-level intersection
    // testing. We loop over the entire sorted set of (key,value) pairs
    // creating Candidate entries for pairs with the same key

    // 1. Soft merge of items (assuming items are sorted)
    size_t num_items = items0.size() + items1.size();
    std::vector<long> merged_item_indices;
    merged_item_indices.reserve(num_items);
    {
        long i = 0, j = 0;
        while (i < items0.size() && j < items1.size()) {
            if (items0[i] < items1[j]) {
                merged_item_indices.push_back(-(i++) - 1);
            } else {
                merged_item_indices.push_back(j++);
            }
        }
        while (i < items0.size()) {
            merged_item_indices.push_back(-(i++) - 1);
        }
        while (j < items1.size()) {
            merged_item_indices.push_back(j++);
        }
    }
    assert(merged_item_indices.size() == num_items);

    const auto get_item = [&](long i) -> const HashItem& {
        return i < 0 ? items0[-(i + 1)] : items1[i];
    };

    // 2. Enumerate hash collisions
#ifdef IPC_TOOLKIT_HASH_GRID_USE_SORT_UNIQUE
    tbb::enumerable_thread_specific<std::vector<Candidate>> storage;
#else
    tbb::enumerable_thread_specific<unordered_set<Candidate>> storage;
#endif

    tbb::parallel_for(
        tbb::blocked_range2d<long>(0l, num_items - 1, 0l, num_items),
        [&](const tbb::blocked_range2d<long>& r) {
            auto& local_candidates = storage.local();

            // i < j
            long i_end = std::min(r.rows().end(), r.cols().end());
            for (long i = r.rows().begin(); i < i_end; i++) {
                const long idx0 = merged_item_indices[i];
                const HashItem& item0 = get_item(idx0);

                // i < r.cols().end() → i + 1 <= r.cols().end()
                long j_begin = std::max(i + 1, r.cols().begin());
                for (long j = j_begin; j < r.cols().end(); j++) {
                    const long idx1 = merged_item_indices[j];
                    const HashItem& item1 = get_item(idx1);

                    if (item0.key != item1.key) {
                        break; // This avoids a brute force comparison
                    }

                    long id0 = item0.id, id1 = item1.id;
                    if (idx0 >= 0 && idx1 < 0) {
                        std::swap(id0, id1);
                    } else if (idx0 >= 0 || idx1 < 0) {
                        continue;
                    }
                    assert(id0 < boxes0.size() && id1 < boxes1.size());

                    if (!can_collide(id0, id1)) {
                        continue;
                    }

                    if (boxes0[id0].intersects(boxes1[id1])) {
#ifdef IPC_TOOLKIT_HASH_GRID_USE_SORT_UNIQUE
                        local_candidates.emplace_back(id0, id1);
#else
                        local_candidates.emplace(id0, id1);
#endif
                    }
                }
            }
        });

#ifdef IPC_TOOLKIT_HASH_GRID_USE_SORT_UNIQUE
    merge_thread_local_vectors(storage, candidates);

    // Remove the duplicate candidates
    tbb::parallel_sort(candidates.begin(), candidates.end());
    auto new_end = std::unique(candidates.begin(), candidates.end());
    candidates.erase(new_end, candidates.end());
#else
    unordered_set<Candidate> candidates_set;
    merge_thread_local_unordered_sets(storage, candidates_set);

    candidates.reserve(candidates_set.size());
    candidates.insert(
        candidates.end(), candidates_set.begin(), candidates_set.end());
#endif
}

template <typename Candidate>
void HashGrid::detect_candidates(
    const std::vector<HashItem>& items,
    const std::vector<AABB>& boxes,
    const std::function<bool(size_t, size_t)>& can_collide,
    std::vector<Candidate>& candidates) const
{
    // Entries with the same key means they share a cell (that cell index
    // hashes to the same key) and should be flagged for low-level
    // intersection testing. So we loop over the entire sorted set of
    // (key,value) pairs creating Candidate entries for pairs with the same key

#ifdef IPC_TOOLKIT_HASH_GRID_USE_SORT_UNIQUE
    tbb::enumerable_thread_specific<std::vector<Candidate>> storage;
#else
    tbb::enumerable_thread_specific<unordered_set<Candidate>> storage;
#endif

    tbb::parallel_for(
        tbb::blocked_range2d<long>(0l, items.size() - 1, 0l, items.size()),
        [&](const tbb::blocked_range2d<long>& r) {
            auto& local_candidates = storage.local();

            // i < j
            long i_end = std::min(r.rows().end(), r.cols().end());
            for (long i = r.rows().begin(); i < i_end; i++) {
                const HashItem& item0 = items[i];
                const AABB& box0 = boxes[item0.id];

                // i < r.cols().end() → i + 1 <= r.cols().end()
                long j_begin = std::max(i + 1, r.cols().begin());
                assert(j_begin > i);
                for (long j = j_begin; j < r.cols().end(); j++) {
                    const HashItem& item1 = items[j];

                    if (item0.key != item1.key) {
                        break; // This avoids a brute force comparison
                    }

                    if (!can_collide(item0.id, item1.id)) {
                        continue;
                    }

                    const AABB& box1 = boxes[item1.id];
                    if (box0.intersects(box1)) {
#ifdef IPC_TOOLKIT_HASH_GRID_USE_SORT_UNIQUE
                        local_candidates.emplace_back(item0.id, item1.id);
#else
                        local_candidates.emplace(item0.id, item1.id);
#endif
                    }
                }
            }
        });

#ifdef IPC_TOOLKIT_HASH_GRID_USE_SORT_UNIQUE
    merge_thread_local_vectors(storage, candidates);

    // Remove the duplicate candidates
    tbb::parallel_sort(candidates.begin(), candidates.end());
    auto new_end = std::unique(candidates.begin(), candidates.end());
    candidates.erase(new_end, candidates.end());
#else
    unordered_set<Candidate> candidates_set;
    merge_thread_local_unordered_sets(storage, candidates_set);

    candidates.reserve(candidates_set.size());
    candidates.insert(
        candidates.end(), candidates_set.begin(), candidates_set.end());
#endif
}

void HashGrid::detect_vertex_vertex_candidates(
    std::vector<VertexVertexCandidate>& candidates) const
{
    detect_candidates(
        vertex_items, vertex_boxes, can_vertices_collide, candidates);
}

void HashGrid::detect_edge_vertex_candidates(
    std::vector<EdgeVertexCandidate>& candidates) const
{
    detect_candidates(
        edge_items, vertex_items, edge_boxes, vertex_boxes,
        [&](size_t ei, size_t vi) { return can_edge_vertex_collide(ei, vi); },
        candidates);
}

void HashGrid::detect_edge_edge_candidates(
    std::vector<EdgeEdgeCandidate>& candidates) const
{
    detect_candidates(
        edge_items, edge_boxes,
        [&](size_t eai, size_t ebi) { return can_edges_collide(eai, ebi); },
        candidates);
}

void HashGrid::detect_face_vertex_candidates(
    std::vector<FaceVertexCandidate>& candidates) const
{
    detect_candidates(
        face_items, vertex_items, face_boxes, vertex_boxes,
        [&](size_t fi, size_t vi) { return can_face_vertex_collide(fi, vi); },
        candidates);
}

void HashGrid::detect_edge_face_candidates(
    std::vector<EdgeFaceCandidate>& candidates) const
{
    detect_candidates(
        edge_items, face_items, edge_boxes, face_boxes,
        [&](size_t ei, size_t fi) { return can_edge_face_collide(ei, fi); },
        candidates);
}

} // namespace ipc
