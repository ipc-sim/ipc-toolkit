#include <ipc/broad_phase/hash_grid.hpp>

#include <algorithm> // std::min/max

#include <tbb/enumerable_thread_specific.h>
#include <tbb/blocked_range2d.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>

#include <ipc/broad_phase/voxel_size_heuristic.hpp>
#include <ipc/utils/vertex_to_min_edge.hpp>
#include <ipc/utils/logger.hpp>

// NOTE: Uncomment for expiermental parallel version of getPairs()
#define IPC_HASH_GRID_PARALLEL_GET_PAIR

namespace ipc {

AABB::AABB(const ArrayMax3d& min, const ArrayMax3d& max)
    : min(min)
    , max(max)
{
    assert(min.size() == max.size());
    assert((min <= max).all());
    // half_extent = (max - min) / 2;
    // center = min + half_extent;
}

bool AABB::are_overlapping(const AABB& a, const AABB& b)
{
    assert(a.min.size() == b.min.size());

    // NOTE: This is a faster check (https://gamedev.stackexchange.com/a/587),
    // but it is inexact and can result in false negatives.
    // return (abs(a.center.x() - b.center.x())
    //         <= (a.half_extent.x() + b.half_extent.x()))
    //     && (abs(a.center.y() - b.center.y())
    //         <= (a.half_extent.y() + b.half_extent.y()))
    //     && (a.min.size() == 2
    //         || abs(a.center.z() - b.center.z())
    //             <= (a.half_extent.z() + b.half_extent.z()));

    // This on the otherhand, is exact because there is no rounding.
    return (a.min <= b.max).all() && (b.min <= a.max).all();
};

typedef tbb::enumerable_thread_specific<std::vector<HashItem>>
    ThreadSpecificHashItems;

void merge_local_items(
    const ThreadSpecificHashItems& storages, std::vector<HashItem>& items)
{
    // size up the hash items
    size_t num_items = items.size();
    for (const auto& local_items : storages) {
        num_items += local_items.size();
    }
    // serial merge!
    items.reserve(num_items);
    for (const auto& local_items : storages) {
        items.insert(items.end(), local_items.begin(), local_items.end());
    }
}

void HashGrid::resizeFromBox(
    const ArrayMax3d& min, const ArrayMax3d& max, double cellSize)
{
    clear();
    assert(cellSize != 0.0);
    m_cellSize = cellSize;
    m_domainMin = min;
    m_domainMax = max;
    m_gridSize = ((max - min) / m_cellSize).ceil().cast<int>().max(1);
    IPC_LOG(trace(
        "hash-grid resized with a size of {:d}x{:d}x{:d}", m_gridSize[0],
        m_gridSize[1], m_gridSize.size() == 3 ? m_gridSize[2] : 1));
}

/// @brief Compute an AABB around a given 2D mesh.
void calculate_mesh_extents(
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    ArrayMax3d& lower_bound,
    ArrayMax3d& upper_bound)
{
    ArrayMax3d lower_bound_t0 = vertices_t0.colwise().minCoeff();
    ArrayMax3d upper_bound_t0 = vertices_t0.colwise().maxCoeff();
    ArrayMax3d lower_bound_t1 = vertices_t1.colwise().minCoeff();
    ArrayMax3d upper_bound_t1 = vertices_t1.colwise().maxCoeff();
    lower_bound = lower_bound_t0.min(lower_bound_t1);
    upper_bound = upper_bound_t0.max(upper_bound_t1);
}

void HashGrid::resize(
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const Eigen::MatrixXi& edges,
    const double inflation_radius)
{
    ArrayMax3d mesh_min, mesh_max;
    calculate_mesh_extents(vertices_t0, vertices_t1, mesh_min, mesh_max);
    double cell_size = suggest_good_voxel_size(
        vertices_t0, vertices_t1, edges, inflation_radius);
    resizeFromBox(
        mesh_min - inflation_radius, mesh_max + inflation_radius, cell_size);
}

/// @brief Compute a AABB for a vertex moving through time (i.e. temporal edge).
void calculate_vertex_extents(
    const VectorMax3d& vertex_t0,
    const VectorMax3d& vertex_t1,
    ArrayMax3d& lower_bound,
    ArrayMax3d& upper_bound)
{
    lower_bound = vertex_t0.cwiseMin(vertex_t1);
    upper_bound = vertex_t0.cwiseMax(vertex_t1);
}

void HashGrid::addVertex(
    const VectorMax3d& vertex_t0,
    const VectorMax3d& vertex_t1,
    const long index,
    const double inflation_radius)
{
    addVertex(vertex_t0, vertex_t1, index, m_vertexItems, inflation_radius);
}

void HashGrid::addVertex(
    const VectorMax3d& vertex_t0,
    const VectorMax3d& vertex_t1,
    const long index,
    std::vector<HashItem>& vertex_items,
    const double inflation_radius) const
{
    ArrayMax3d lower_bound, upper_bound;
    calculate_vertex_extents(vertex_t0, vertex_t1, lower_bound, upper_bound);
    addElement(
        AABB(lower_bound - inflation_radius, upper_bound + inflation_radius),
        index, vertex_items);
}

void HashGrid::addVertices(
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const double inflation_radius)
{
    assert(vertices_t0.rows() == vertices_t1.rows());

    ThreadSpecificHashItems storage;

    tbb::parallel_for(
        tbb::blocked_range<long>(0l, long(vertices_t0.rows())),
        [&](const tbb::blocked_range<long>& range) {
            ThreadSpecificHashItems::reference local_items = storage.local();

            for (long i = range.begin(); i != range.end(); i++) {
                addVertex(
                    vertices_t0.row(i), vertices_t1.row(i), i, local_items,
                    inflation_radius);
            }
        });

    merge_local_items(storage, m_vertexItems);
}

void HashGrid::addSelectVertices(
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const Eigen::VectorXi& vertex_indices,
    const double inflation_radius)
{
    assert(vertices_t0.rows() == vertices_t1.rows());

    ThreadSpecificHashItems storage;

    tbb::parallel_for(
        tbb::blocked_range<long>(0l, long(vertex_indices.size())),
        [&](const tbb::blocked_range<long>& range) {
            ThreadSpecificHashItems::reference local_items = storage.local();

            for (long i = range.begin(); i != range.end(); i++) {
                const size_t vi = vertex_indices[i];
                addVertex(
                    vertices_t0.row(vi), vertices_t1.row(vi), vi, local_items,
                    inflation_radius);
            }
        });

    merge_local_items(storage, m_vertexItems);
}

void HashGrid::addVerticesFromEdges(
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const Eigen::MatrixXi& edges,
    const double inflation_radius)
{
    assert(vertices_t0.rows() == vertices_t1.rows());

    std::vector<size_t> V2E = vertex_to_min_edge(vertices_t0.rows(), edges);

    ThreadSpecificHashItems storage;

    tbb::parallel_for(
        tbb::blocked_range<long>(0l, long(edges.rows())),
        [&](const tbb::blocked_range<long>& range) {
            ThreadSpecificHashItems::reference local_items = storage.local();

            for (long ei = range.begin(); ei < range.end(); ei++) {
                for (long ej = 0; ej < edges.cols(); ej++) {
                    const size_t vi = edges(ei, ej);
                    if (V2E[vi] == ei) {
                        addVertex(
                            vertices_t0.row(vi), vertices_t1.row(vi), vi,
                            local_items, inflation_radius);
                    }
                }
            }
        });

    merge_local_items(storage, m_vertexItems);
}

/// @brief Compute a AABB for an edge moving through time (i.e. temporal quad).
void calculate_edge_extents(
    const VectorMax3d& edge_vertex0_t0,
    const VectorMax3d& edge_vertex1_t0,
    const VectorMax3d& edge_vertex0_t1,
    const VectorMax3d& edge_vertex1_t1,
    ArrayMax3d& lower_bound,
    ArrayMax3d& upper_bound)
{
    lower_bound = edge_vertex0_t0.cwiseMin(edge_vertex1_t0)
                      .cwiseMin(edge_vertex0_t1)
                      .cwiseMin(edge_vertex1_t1);
    upper_bound = edge_vertex0_t0.cwiseMax(edge_vertex1_t0)
                      .cwiseMax(edge_vertex0_t1)
                      .cwiseMax(edge_vertex1_t1);
}

void HashGrid::addEdge(
    const VectorMax3d& edge_vertex0_t0,
    const VectorMax3d& edge_vertex1_t0,
    const VectorMax3d& edge_vertex0_t1,
    const VectorMax3d& edge_vertex1_t1,
    const long index,
    const double inflation_radius)
{
    addEdge(
        edge_vertex0_t0, edge_vertex1_t0, edge_vertex0_t1, edge_vertex1_t1,
        index, m_edgeItems, inflation_radius);
}

void HashGrid::addEdge(
    const VectorMax3d& edge_vertex0_t0,
    const VectorMax3d& edge_vertex1_t0,
    const VectorMax3d& edge_vertex0_t1,
    const VectorMax3d& edge_vertex1_t1,
    const long index,
    std::vector<HashItem>& edge_items,
    const double inflation_radius) const
{
    ArrayMax3d lower_bound, upper_bound;
    calculate_edge_extents(
        edge_vertex0_t0, edge_vertex1_t0, edge_vertex0_t1, edge_vertex1_t1,
        lower_bound, upper_bound);
    addElement(
        AABB(lower_bound - inflation_radius, upper_bound + inflation_radius),
        index, edge_items);
}

void HashGrid::addEdges(
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const Eigen::MatrixXi& edges,
    const double inflation_radius)
{
    assert(vertices_t0.rows() == vertices_t1.rows());

    ThreadSpecificHashItems storage;

    tbb::parallel_for(
        tbb::blocked_range<long>(0l, long(edges.rows())),
        [&](const tbb::blocked_range<long>& range) {
            ThreadSpecificHashItems::reference local_items = storage.local();

            for (long i = range.begin(); i != range.end(); i++) {
                addEdge(
                    vertices_t0.row(edges(i, 0)), vertices_t0.row(edges(i, 1)),
                    vertices_t1.row(edges(i, 0)), vertices_t1.row(edges(i, 1)),
                    i, local_items, inflation_radius);
            }
        });

    merge_local_items(storage, m_edgeItems);
}

/// @brief Compute a AABB for an edge moving through time (i.e. temporal quad).
void calculate_face_extents(
    const VectorMax3d& face_vertex0_t0,
    const VectorMax3d& face_vertex1_t0,
    const VectorMax3d& face_vertex2_t0,
    const VectorMax3d& face_vertex0_t1,
    const VectorMax3d& face_vertex1_t1,
    const VectorMax3d& face_vertex2_t1,
    ArrayMax3d& lower_bound,
    ArrayMax3d& upper_bound)
{
    lower_bound = face_vertex0_t0.cwiseMin(face_vertex1_t0)
                      .cwiseMin(face_vertex2_t0)
                      .cwiseMin(face_vertex0_t1)
                      .cwiseMin(face_vertex1_t1)
                      .cwiseMin(face_vertex2_t1);
    upper_bound = face_vertex0_t0.cwiseMax(face_vertex1_t0)
                      .cwiseMax(face_vertex2_t0)
                      .cwiseMax(face_vertex0_t1)
                      .cwiseMax(face_vertex1_t1)
                      .cwiseMax(face_vertex2_t1);
}

void HashGrid::addFace(
    const VectorMax3d& face_vertex0_t0,
    const VectorMax3d& face_vertex1_t0,
    const VectorMax3d& face_vertex2_t0,
    const VectorMax3d& face_vertex0_t1,
    const VectorMax3d& face_vertex1_t1,
    const VectorMax3d& face_vertex2_t1,
    const long index,
    const double inflation_radius)
{
    addFace(
        face_vertex0_t0, face_vertex1_t0, face_vertex2_t0, face_vertex0_t1,
        face_vertex1_t1, face_vertex2_t1, index, m_faceItems, inflation_radius);
}

void HashGrid::addFace(
    const VectorMax3d& face_vertex0_t0,
    const VectorMax3d& face_vertex1_t0,
    const VectorMax3d& face_vertex2_t0,
    const VectorMax3d& face_vertex0_t1,
    const VectorMax3d& face_vertex1_t1,
    const VectorMax3d& face_vertex2_t1,
    const long index,
    std::vector<HashItem>& face_items,
    const double inflation_radius) const
{
    ArrayMax3d lower_bound, upper_bound;
    calculate_face_extents(
        face_vertex0_t0, face_vertex1_t0, face_vertex2_t0, //
        face_vertex0_t1, face_vertex1_t1, face_vertex2_t1, //
        lower_bound, upper_bound);
    addElement(
        AABB(lower_bound - inflation_radius, upper_bound + inflation_radius),
        index, face_items);
}

void HashGrid::addFaces(
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const Eigen::MatrixXi& faces,
    const double inflation_radius)
{
    assert(vertices_t0.rows() == vertices_t1.rows());

    ThreadSpecificHashItems storage;

    tbb::parallel_for(
        tbb::blocked_range<long>(0l, long(faces.rows())),
        [&](const tbb::blocked_range<long>& range) {
            ThreadSpecificHashItems::reference local_items = storage.local();

            for (long i = range.begin(); i != range.end(); i++) {
                addFace(
                    vertices_t0.row(faces(i, 0)), vertices_t0.row(faces(i, 1)),
                    vertices_t0.row(faces(i, 2)), vertices_t1.row(faces(i, 0)),
                    vertices_t1.row(faces(i, 1)), vertices_t1.row(faces(i, 2)),
                    i, local_items, inflation_radius);
            }
        });

    merge_local_items(storage, m_faceItems);
}

void HashGrid::addElement(
    const AABB& aabb, const int id, std::vector<HashItem>& items) const
{
    ArrayMax3i int_min =
        ((aabb.getMin() - m_domainMin) / m_cellSize).cast<int>();
    // We can round down to -1, but not less
    assert((int_min >= -1).all());
    assert((int_min <= m_gridSize).all());
    int_min = int_min.max(0).min(m_gridSize - 1);

    ArrayMax3i int_max =
        ((aabb.getMax() - m_domainMin) / m_cellSize).cast<int>();
    assert((int_max >= -1).all());
    assert((int_max <= m_gridSize).all());
    int_max = int_max.max(0).min(m_gridSize - 1);
    assert((int_min <= int_max).all());

    int min_z = int_min.size() == 3 ? int_min.z() : 0;
    int max_z = int_max.size() == 3 ? int_max.z() : 0;
    for (int x = int_min.x(); x <= int_max.x(); ++x) {
        for (int y = int_min.y(); y <= int_max.y(); ++y) {
            for (int z = min_z; z <= max_z; ++z) {
                items.emplace_back(hash(x, y, z), id, aabb);
            }
        }
    }
}

template <typename T>
void HashGrid::getPairs(
    const std::function<bool(int, int)>& is_endpoint,
    const std::function<bool(int, int)>& can_collide,
    std::vector<HashItem>& items0,
    std::vector<HashItem>& items1,
    T& candidates)
{
    // Sorted all they (key,value) pairs, where key is the hash key, and value
    // is the element index
    tbb::parallel_sort(items0.begin(), items0.end());
    tbb::parallel_sort(items1.begin(), items1.end());

    // Entries with the same key means they share a cell (that cell index
    // hashes to the same key) and should be flagged for low-level intersection
    // testing. We loop over the entire sorted set of (key,value) pairs
    // creating Candidate entries for pairs with the same key

    // 1. Find start of cells in items0 and items1. Store the start and end of
    // each cell in a vector.
    int num_cells = m_gridSize.prod();

    std::vector<long> cell_starts0(num_cells, -1);
    std::vector<long> cell_ends0(num_cells, -1);
    for (long i = 0; i < items0.size(); i++) {
        if (i == 0 || items0[i].key != items0[i - 1].key) {
            cell_starts0[items0[i].key] = i;
        }
        if (i == items0.size() - 1 || items0[i].key != items0[i + 1].key) {
            cell_ends0[items0[i].key] = i;
        }
    }

    std::vector<long> cell_starts1(num_cells, -1);
    std::vector<long> cell_ends1(num_cells, -1);
    for (long i = 0; i < items1.size(); i++) {
        if (i == 0 || items1[i].key != items1[i - 1].key) {
            cell_starts1[items1[i].key] = i;
        }
        if (i == items1.size() - 1 || items1[i].key != items1[i + 1].key) {
            cell_ends1[items1[i].key] = i;
        }
    }

    // 2. Perform a double for loop for each cell, building the candidates.
#ifdef IPC_HASH_GRID_PARALLEL_GET_PAIR
    tbb::enumerable_thread_specific<T> storages;
    tbb::parallel_for(
        tbb::blocked_range<int>(0, num_cells),
        [&](const tbb::blocked_range<int>& c_ids) {
            for (int ci = c_ids.begin(); ci < c_ids.end(); ci++) {
                auto& local_candidates = storages.local();
#else
    for (int ci = 0; ci < num_cells; ci++) {
        auto& local_candidates = candidates;
#endif
                // Check if the cell is empty
                if (cell_starts0[ci] < 0 || cell_starts1[ci] < 0) {
                    continue;
                }
                assert(cell_ends0[ci] >= 0 && cell_ends1[ci] >= 0);

                for (int i = cell_starts0[ci]; i <= cell_ends0[ci]; i++) {
                    for (int j = cell_starts1[ci]; j <= cell_ends1[ci]; j++) {
                        if (!is_endpoint(items0[i].id, items1[j].id)
                            && can_collide(items0[i].id, items1[j].id)
                            && AABB::are_overlapping(
                                items0[i].aabb, items1[j].aabb)) {
                            local_candidates.emplace_back(
                                items0[i].id, items1[j].id);
                        }
                    }
                }
            }
#ifdef IPC_HASH_GRID_PARALLEL_GET_PAIR
        });

    // size up the candidates
    size_t num_candidates = candidates.size();
    for (const auto& local_candidates : storages) {
        num_candidates += local_candidates.size();
    }
    // serial merge!
    candidates.reserve(num_candidates);
    for (const auto& local_candidates : storages) {
        candidates.insert(
            candidates.end(), local_candidates.begin(), local_candidates.end());
    }
#endif

    // WARNING: This version is buggy and I am not sure why.
    /*
    int i = 0, j_start = 0;
    while (i < items0.size() && j_start < items1.size()) {
        const HashItem& item0 = items0[i];

        int j = j_start;
        while (j < items1.size()) {
            const HashItem& item1 = items1[j];

            if (item0.key == item1.key) {
                if (!is_endpoint(item0.id, item1.id)
                    && can_collide(item0.id, item1.id)
                    && AABB::are_overlapping(item0.aabb, item1.aabb)) {
                    candidates.emplace_back(item0.id, item1.id);
                }
            } else {
                break;
            }
            j++;
        }

        if (i == items0.size() - 1 || item0.key != items0[i + 1].key) {
            j_start = j;
        }
        i++;
    }
    */

    // Remove the duplicate candidates
    tbb::parallel_sort(candidates.begin(), candidates.end());
    auto new_end = std::unique(candidates.begin(), candidates.end());
    candidates.erase(new_end, candidates.end());
}

template <typename T>
void HashGrid::getPairs(
    const std::function<bool(int, int)>& is_endpoint,
    const std::function<bool(int, int)>& can_collide,
    std::vector<HashItem>& items,
    T& candidates)
{
    // Sorted all they (key,value) pairs, where key is the hash key, and
    // value is the element index
    tbb::parallel_sort(items.begin(), items.end());

#ifndef IPC_HASH_GRID_PARALLEL_GET_PAIR
    // Entries with the same key means they share a cell (that cell index
    // hashes to the same key) and should be flagged for low-level
    // intersection testing. So we loop over the entire sorted set of
    // (key,value) pairs creating Candidate entries for pairs with the same key
    for (int i = 0; i < items.size(); i++) {
        const HashItem& item0 = items[i];
        for (int j = i + 1; j < items.size(); j++) {
            const HashItem& item1 = items[j];
            if (item0.key == item1.key) {
                if (!is_endpoint(item0.id, item1.id)
                    && can_collide(item0.id, item1.id)
                    && AABB::are_overlapping(item0.aabb, item1.aabb)) {
                    candidates.emplace_back(item0.id, item1.id);
                }
            } else {
                break; // This avoids a brute force comparison
            }
        }
    }
#else
    tbb::enumerable_thread_specific<T> storages;

    tbb::parallel_for(
        tbb::blocked_range2d<int>(0, items.size(), 0, items.size()),
        [&](const tbb::blocked_range2d<int>& r) {
            if (r.rows().begin() >= r.cols().end()) {
                return; // i needs to be less than j
            }

            auto& local_candidates = storages.local();

            for (int i = r.rows().begin(); i < r.rows().end(); i++) {
                const HashItem& item0 = items[i];

                if (i >= r.cols().end()) {
                    return; // i will increase but r.cols().end() will not
                }

                // i < r.cols().end() → i + 1 <= r.cols().end()
                int j_start = std::max(i + 1, r.cols().begin());
                assert(j_start > i);

                for (int j = j_start; j < r.cols().end(); j++) {
                    const HashItem& item1 = items[j];

                    if (item0.key == item1.key) {
                        if (!is_endpoint(item0.id, item1.id)
                            && can_collide(item0.id, item1.id)
                            && AABB::are_overlapping(item0.aabb, item1.aabb)) {
                            local_candidates.emplace_back(item0.id, item1.id);
                        }
                    } else {
                        break; // This avoids a brute force comparison
                    }
                }
            }
        });

    // size up the candidates
    size_t num_candidates = candidates.size();
    for (const auto& local_candidates : storages) {
        num_candidates += local_candidates.size();
    }
    // serial merge!
    candidates.reserve(num_candidates);
    for (const auto& local_candidates : storages) {
        candidates.insert(
            candidates.end(), local_candidates.begin(), local_candidates.end());
    }
#endif

    // Remove the duplicate candidates
    tbb::parallel_sort(candidates.begin(), candidates.end());
    auto new_end = std::unique(candidates.begin(), candidates.end());
    candidates.erase(new_end, candidates.end());
}

void HashGrid::getVertexEdgePairs(
    const Eigen::MatrixXi& edges,
    std::vector<EdgeVertexCandidate>& ev_candidates,
    const std::function<bool(size_t, size_t)>& can_vertices_collide)
{
    auto is_endpoint = [&](int ei, int vi) {
        return edges(ei, 0) == vi || edges(ei, 1) == vi;
    };

    auto can_collide = [&](int ei, int vi) {
        return can_vertices_collide(vi, edges(ei, 0))
            || can_vertices_collide(vi, edges(ei, 1));
    };

    getPairs(
        is_endpoint, can_collide, m_edgeItems, m_vertexItems, ev_candidates);
}

void HashGrid::getEdgeEdgePairs(
    const Eigen::MatrixXi& edges,
    std::vector<EdgeEdgeCandidate>& ee_candidates,
    const std::function<bool(size_t, size_t)>& can_vertices_collide)
{
    auto is_endpoint = [&](int ei, int ej) {
        return edges(ei, 0) == edges(ej, 0) || edges(ei, 0) == edges(ej, 1)
            || edges(ei, 1) == edges(ej, 0) || edges(ei, 1) == edges(ej, 1);
    };

    auto can_collide = [&](int ei, int ej) {
        return can_vertices_collide(edges(ei, 0), edges(ej, 0))
            || can_vertices_collide(edges(ei, 0), edges(ej, 1))
            || can_vertices_collide(edges(ei, 1), edges(ej, 0))
            || can_vertices_collide(edges(ei, 1), edges(ej, 1));
    };

    getPairs(is_endpoint, can_collide, m_edgeItems, ee_candidates);
}

void HashGrid::getEdgeFacePairs(
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    std::vector<EdgeFaceCandidate>& ef_candidates,
    const std::function<bool(size_t, size_t)>& can_vertices_collide)
{
    auto is_endpoint = [&](int ei, int fi) {
        // Check if the edge and face have a common end-point
        return edges(ei, 0) == faces(fi, 0) || edges(ei, 0) == faces(fi, 1)
            || edges(ei, 0) == faces(fi, 2) || edges(ei, 1) == faces(fi, 0)
            || edges(ei, 1) == faces(fi, 1) || edges(ei, 1) == faces(fi, 2);
    };

    auto can_collide = [&](int ei, int fi) {
        return can_vertices_collide(edges(ei, 0), faces(fi, 0))
            || can_vertices_collide(edges(ei, 0), faces(fi, 1))
            || can_vertices_collide(edges(ei, 0), faces(fi, 2))
            || can_vertices_collide(edges(ei, 1), faces(fi, 0))
            || can_vertices_collide(edges(ei, 1), faces(fi, 1))
            || can_vertices_collide(edges(ei, 1), faces(fi, 2));
    };

    getPairs(is_endpoint, can_collide, m_edgeItems, m_faceItems, ef_candidates);
}

void HashGrid::getFaceVertexPairs(
    const Eigen::MatrixXi& faces,
    std::vector<FaceVertexCandidate>& fv_candidates,
    const std::function<bool(size_t, size_t)>& can_vertices_collide)
{
    auto is_endpoint = [&](int fi, int vi) {
        return vi == faces(fi, 0) || vi == faces(fi, 1) || vi == faces(fi, 2);
    };

    auto can_collide = [&](int fi, int vi) {
        return can_vertices_collide(vi, faces(fi, 0))
            || can_vertices_collide(vi, faces(fi, 1))
            || can_vertices_collide(vi, faces(fi, 2));
    };

    getPairs(
        is_endpoint, can_collide, m_faceItems, m_vertexItems, fv_candidates);
}

} // namespace ipc
