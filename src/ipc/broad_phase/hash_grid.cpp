#include "hash_grid.hpp"

#include <ipc/broad_phase/voxel_size_heuristic.hpp>
#include <ipc/utils/logger.hpp>
#include <ipc/utils/merge_thread_local.hpp>

#include <tbb/blocked_range2d.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>

#include <algorithm> // std::min/max

#define IPC_TOOLKIT_HASH_GRID_USE_SORT_UNIQUE // else use unordered_set

using namespace std::placeholders;

namespace ipc {

void HashGrid::build(
    Eigen::ConstRef<Eigen::MatrixXi> edges,
    Eigen::ConstRef<Eigen::MatrixXi> faces)
{
    BroadPhase::build(edges, faces);

    Eigen::Array3d mesh_min, mesh_max;
    compute_mesh_aabb(mesh_min, mesh_max);

    const double cell_size =
        suggest_good_voxel_size(edges.rows() > 0 ? edge_boxes : vertex_boxes);
    resize(mesh_min, mesh_max, cell_size);

    insert_boxes();
}

void HashGrid::resize(
    Eigen::ConstRef<Eigen::Array3d> domain_min,
    Eigen::ConstRef<Eigen::Array3d> domain_max,
    const double cell_size)
{
    assert(cell_size > 0.0);
    assert(std::isfinite(cell_size));

    m_domain_min = domain_min;
    m_domain_max = domain_max;
    m_cell_size = cell_size;
    m_grid_size =
        ((domain_max - domain_min) / cell_size).ceil().cast<int>().max(1);

    logger().trace(
        "hash-grid resized with a size of {:d}x{:d}x{:d}", grid_size()[0],
        grid_size()[1], grid_size()[2]);
}

void HashGrid::insert_boxes()
{
    insert_boxes(this->vertex_boxes, vertex_items);
    insert_boxes(this->edge_boxes, edge_items);
    insert_boxes(this->face_boxes, face_items);
}

void HashGrid::insert_boxes(
    const AABBs& boxes, std::vector<HashItem>& items) const
{
    tbb::enumerable_thread_specific<std::vector<HashItem>> storage;

    tbb::parallel_for(0L, long(boxes.size()), [&](long i) {
        insert_box(boxes[i], i, storage.local());
    });

    merge_thread_local_vectors(storage, items);

    // Sorted all they (key, value) pairs, where key is the hash key, and
    // value is the element index
    tbb::parallel_sort(items.begin(), items.end());
}

void HashGrid::insert_box(
    const AABB& aabb, const long id, std::vector<HashItem>& items) const
{
    Eigen::Array3i int_min =
        ((aabb.min - domain_min()) / cell_size()).cast<int>();
    // We can round down to -1, but not less
    assert((int_min >= -1).all());
    assert((int_min <= grid_size()).all());
    int_min = int_min.max(0).min(grid_size() - 1);

    Eigen::Array3i int_max =
        ((aabb.max - domain_min()) / cell_size()).cast<int>();
    assert((int_max >= -1).all());
    assert((int_max <= grid_size()).all());
    int_max = int_max.max(0).min(grid_size() - 1);
    assert((int_min <= int_max).all());

    const int min_z = dim == 3 ? int_min.z() : 0;
    const int max_z = dim == 3 ? int_max.z() : 0;
    for (int x = int_min.x(); x <= int_max.x(); ++x) {
        for (int y = int_min.y(); y <= int_max.y(); ++y) {
            for (int z = min_z; z <= max_z; ++z) {
                items.emplace_back(hash(x, y, z), id);
            }
        }
    }
}

namespace {
    // Join sorted cell runs, then enumerate only pairs inside a common cell.
    // Large cells retain a 2D parallel split; small cells share an outer task.
    template <typename Candidate, bool SameItems>
    void detect_cell_pairs(
        const std::vector<HashItem>& items0,
        const std::vector<HashItem>& items1,
        const AABBs& boxes0,
        const AABBs& boxes1,
        const std::function<bool(size_t, size_t)>& can_collide,
        std::vector<Candidate>& candidates)
    {
        struct Cell {
            long begin0, end0, begin1, end1;
        };
        std::vector<Cell> cells;
        const auto run_end = [](const auto& items, long begin) {
            long end = begin + 1;
            while (end < items.size() && items[end].key == items[begin].key)
                ++end;
            return end;
        };
        long i = 0, j = 0;
        while (i < items0.size() && j < items1.size()) {
            if constexpr (SameItems) {
                const long ie = run_end(items0, i);
                if (ie - i > 1)
                    cells.push_back({ i, ie, i, ie });
                i = j = ie;
            } else {
                if (items0[i].key < items1[j].key) {
                    i = run_end(items0, i);
                    continue;
                }
                if (items1[j].key < items0[i].key) {
                    j = run_end(items1, j);
                    continue;
                }
                const long ie = run_end(items0, i);
                const long je = run_end(items1, j);
                cells.push_back({ i, ie, j, je });
                i = ie;
                j = je;
            }
        }

#ifdef IPC_TOOLKIT_HASH_GRID_USE_SORT_UNIQUE
        tbb::enumerable_thread_specific<std::vector<Candidate>> storage;
#else
        tbb::enumerable_thread_specific<unordered_set<Candidate>> storage;
#endif
        const auto emit = [&](long begin0, long end0, long begin1, long end1) {
            auto& local_candidates = storage.local();
            if constexpr (SameItems)
                end0 = std::min(end0, end1);
            for (long a = begin0; a < end0; ++a) {
                const auto id0 = items0[a].id;
                const AABB& box0 = boxes0[id0];
                const long first = SameItems ? std::max(a + 1, begin1) : begin1;
                for (long b = first; b < end1; ++b) {
                    const auto id1 = items1[b].id;
                    if (!box0.intersects(boxes1[id1]) || !can_collide(id0, id1))
                        continue;
#ifdef IPC_TOOLKIT_HASH_GRID_USE_SORT_UNIQUE
                    local_candidates.emplace_back(id0, id1);
#else
                    local_candidates.emplace(id0, id1);
#endif
                }
            }
        };
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, cells.size()),
            [&](const auto& range) {
                for (size_t k = range.begin(); k < range.end(); ++k) {
                    const auto& c = cells[k];
                    if (c.end0 - c.begin0 <= 64 && c.end1 - c.begin1 <= 64) {
                        emit(c.begin0, c.end0, c.begin1, c.end1);
                    } else {
                        tbb::parallel_for(
                            tbb::blocked_range2d<long>(
                                c.begin0, c.end0, 64, c.begin1, c.end1, 64),
                            [&](const auto& r) {
                                emit(
                                    r.rows().begin(), r.rows().end(),
                                    r.cols().begin(), r.cols().end());
                            });
                    }
                }
            });
#ifdef IPC_TOOLKIT_HASH_GRID_USE_SORT_UNIQUE
        merge_thread_local_vectors(storage, candidates);
        tbb::parallel_sort(candidates.begin(), candidates.end());
        candidates.erase(
            std::unique(candidates.begin(), candidates.end()),
            candidates.end());
#else
        unordered_set<Candidate> candidates_set;
        merge_thread_local_unordered_sets(storage, candidates_set);
        candidates.reserve(candidates_set.size());
        candidates.insert(
            candidates.end(), candidates_set.begin(), candidates_set.end());
#endif
    }
} // namespace

template <typename Candidate>
void HashGrid::detect_candidates(
    const std::vector<HashItem>& items0,
    const std::vector<HashItem>& items1,
    const AABBs& boxes0,
    const AABBs& boxes1,
    const std::function<bool(size_t, size_t)>& can_collide,
    std::vector<Candidate>& candidates) const
{
    detect_cell_pairs<Candidate, false>(
        items0, items1, boxes0, boxes1, can_collide, candidates);
}

template <typename Candidate>
void HashGrid::detect_candidates(
    const std::vector<HashItem>& items,
    const AABBs& boxes,
    const std::function<bool(size_t, size_t)>& can_collide,
    std::vector<Candidate>& candidates) const
{
    detect_cell_pairs<Candidate, true>(
        items, items, boxes, boxes, can_collide, candidates);
}

void HashGrid::detect_vertex_vertex_candidates(
    std::vector<VertexVertexCandidate>& candidates) const
{
    candidates.clear();
    detect_candidates(
        vertex_items, vertex_boxes, can_vertices_collide, candidates);
}

void HashGrid::detect_edge_vertex_candidates(
    std::vector<EdgeVertexCandidate>& candidates) const
{
    candidates.clear();
    detect_candidates(
        edge_items, vertex_items, edge_boxes, vertex_boxes,
        std::bind(&HashGrid::can_edge_vertex_collide, this, _1, _2),
        candidates);
}

void HashGrid::detect_edge_edge_candidates(
    std::vector<EdgeEdgeCandidate>& candidates) const
{
    candidates.clear();
    detect_candidates(
        edge_items, edge_boxes,
        std::bind(&HashGrid::can_edges_collide, this, _1, _2), candidates);
}

void HashGrid::detect_face_vertex_candidates(
    std::vector<FaceVertexCandidate>& candidates) const
{
    candidates.clear();
    detect_candidates(
        face_items, vertex_items, face_boxes, vertex_boxes,
        std::bind(&HashGrid::can_face_vertex_collide, this, _1, _2),
        candidates);
}

void HashGrid::detect_edge_face_candidates(
    std::vector<EdgeFaceCandidate>& candidates) const
{
    candidates.clear();
    detect_candidates(
        edge_items, face_items, edge_boxes, face_boxes,
        std::bind(&HashGrid::can_edge_face_collide, this, _1, _2), candidates);
}

void HashGrid::detect_face_face_candidates(
    std::vector<FaceFaceCandidate>& candidates) const
{
    candidates.clear();
    detect_candidates(
        face_items, face_boxes,
        std::bind(&HashGrid::can_faces_collide, this, _1, _2), candidates);
}

} // namespace ipc
