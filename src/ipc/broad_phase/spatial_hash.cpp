// Modified version of SpatialHash.hpp from IPC codebase.
// Originally created by Minchen Li.
#include "spatial_hash.hpp"

#include <ipc/ccd/aabb.hpp>
#include <ipc/broad_phase/voxel_size_heuristic.hpp>
#include <ipc/utils/merge_thread_local.hpp>

#include <ipc/config.hpp>

#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>

using namespace std::placeholders;

namespace ipc {

namespace {
    inline void tbb_parallel_block_range_for(
        const size_t start_i,
        const size_t end_i,
        const std::function<void(size_t)>& body)
    {
        tbb::parallel_for(
            tbb::blocked_range(start_i, end_i),
            [&](const tbb::blocked_range<size_t>& range) {
                for (size_t i = range.begin(); i != range.end(); i++) {
                    body(i);
                }
            });
    }

    void fill_primitive_to_voxels(
        const Eigen::Array3i& min_voxel,
        const Eigen::Array3i& max_voxel,
        const ArrayMax3i& voxel_count,
        const int voxel_count_0x1,
        std::vector<int>& primitive_to_voxels)
    {
        assert((min_voxel <= max_voxel).all());
        assert(voxel_count_0x1 == voxel_count[0] * voxel_count[1]);

        primitive_to_voxels.reserve((max_voxel - min_voxel + 1).prod());

        for (int iz = min_voxel[2]; iz <= max_voxel[2]; iz++) {
            int z_offset = iz * voxel_count_0x1;

            for (int iy = min_voxel[1]; iy <= max_voxel[1]; iy++) {
                int yz_offset = iy * voxel_count[0] + z_offset;

                for (int ix = min_voxel[0]; ix <= max_voxel[0]; ix++) {
                    primitive_to_voxels.emplace_back(ix + yz_offset);
                }
            }
        }
    }

    void fill_voxel_to_primitives(
        const size_t num_primitives,
        const size_t primitive_offset,
        const std::vector<std::vector<int>>& primitive_to_voxels,
        unordered_map<int, std::vector<int>>& voxel_to_primitives)
    {
        for (int i = 0; i < num_primitives; i++) {
            for (const auto& voxel : primitive_to_voxels[i]) {
                voxel_to_primitives[voxel].emplace_back(i + primitive_offset);
            }
        }
    }
} // namespace

void SpatialHash::build(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    double inflation_radius,
    double voxel_size)
{
    build(vertices, vertices, edges, faces, inflation_radius, voxel_size);
}

void SpatialHash::build(
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    double inflation_radius,
    double voxel_size)
{
    const size_t num_vertices = vertices_t0.rows();
    dim = vertices_t0.cols();

    assert(vertices_t1.rows() == num_vertices && vertices_t1.cols() == dim);

    // also calls clear()
    BroadPhase::build(vertices_t0, vertices_t1, edges, faces, inflation_radius);

    built_in_radius = inflation_radius;

    if (voxel_size <= 0) {
        voxel_size = suggest_good_voxel_size(
            vertices_t0, vertices_t1, edges, inflation_radius);
    }

    left_bottom_corner = vertices_t0.colwise().minCoeff().cwiseMin(
        vertices_t1.colwise().minCoeff());
    right_top_corner = vertices_t0.colwise().maxCoeff().cwiseMax(
        vertices_t1.colwise().maxCoeff());

    AABB::conservative_inflation(
        left_bottom_corner, right_top_corner, inflation_radius);

    one_div_voxelSize = 1.0 / voxel_size;

    const ArrayMax3d range = right_top_corner - left_bottom_corner;
    voxel_count = (range * one_div_voxelSize).ceil().template cast<int>();
    if (voxel_count.minCoeff() <= 0) {
        // cast overflow due to huge search direction
        one_div_voxelSize = 1.0 / (range.maxCoeff() * 1.01);
        voxel_count.setOnes();
    }
    voxel_count_0x1 = voxel_count[0] * voxel_count[1];

    // ------------------------------------------------------------------------
    // precompute vertex min and max voxel axis indices

    std::vector<Eigen::Array3i> vertex_min_voxel_axis_index(
        num_vertices, Eigen::Array3i::Zero());
    std::vector<Eigen::Array3i> vertex_max_voxel_axis_index(
        num_vertices, Eigen::Array3i::Zero());
    tbb::parallel_for(size_t(0), num_vertices, [&](size_t vi) {
        ArrayMax3d v_min = vertices_t0.row(vi).cwiseMin(vertices_t1.row(vi));
        ArrayMax3d v_max = vertices_t0.row(vi).cwiseMax(vertices_t1.row(vi));
        AABB::conservative_inflation(v_min, v_max, inflation_radius);

        vertex_min_voxel_axis_index[vi].head(dim) =
            locate_voxel_axis_index(v_min);
        vertex_max_voxel_axis_index[vi].head(dim) =
            locate_voxel_axis_index(v_max);
    });

    // ------------------------------------------------------------------------

    point_to_voxels.resize(num_vertices);
    tbb_parallel_block_range_for(0, num_vertices, [&](size_t vi) {
        fill_primitive_to_voxels(
            vertex_min_voxel_axis_index[vi], vertex_max_voxel_axis_index[vi],
            voxel_count, voxel_count_0x1, point_to_voxels[vi]);
    });

    edge_to_voxels.resize(edges.rows());
    tbb_parallel_block_range_for(0, edges.rows(), [&](size_t ei) {
        fill_primitive_to_voxels(
            vertex_min_voxel_axis_index[edges(ei, 0)].min(
                vertex_min_voxel_axis_index[edges(ei, 1)]),
            vertex_max_voxel_axis_index[edges(ei, 0)].max(
                vertex_max_voxel_axis_index[edges(ei, 1)]),
            voxel_count, voxel_count_0x1, edge_to_voxels[ei]);
    });

    face_to_voxels.resize(faces.rows());
    tbb_parallel_block_range_for(0, faces.rows(), [&](size_t fi) {
        fill_primitive_to_voxels(
            vertex_min_voxel_axis_index[faces(fi, 0)]
                .min(vertex_min_voxel_axis_index[faces(fi, 1)])
                .min(vertex_min_voxel_axis_index[faces(fi, 2)]),
            vertex_max_voxel_axis_index[faces(fi, 0)]
                .max(vertex_max_voxel_axis_index[faces(fi, 1)])
                .max(vertex_max_voxel_axis_index[faces(fi, 2)]),
            voxel_count, voxel_count_0x1, face_to_voxels[fi]);
    });

    // ------------------------------------------------------------------------

    edge_start_ind = num_vertices;
    tri_start_ind = edge_start_ind + edges.rows();

    fill_voxel_to_primitives(
        num_vertices, 0, point_to_voxels, voxel_to_primitives);
    fill_voxel_to_primitives(
        edges.rows(), edge_start_ind, edge_to_voxels, voxel_to_primitives);
    fill_voxel_to_primitives(
        faces.rows(), tri_start_ind, face_to_voxels, voxel_to_primitives);
}

void SpatialHash::query_point_for_points(
    int vi, unordered_set<int>& vert_ids) const
{
    vert_ids.clear();
    for (const int voxel : point_to_voxels[vi]) {
        for (const int id : voxel_to_primitives.at(voxel)) {
            if (is_vertex_index(id) && id > vi) {
                vert_ids.insert(id);
            }
        }
    }
}

void SpatialHash::query_point_for_edges(
    int vi, unordered_set<int>& edge_ids) const
{
    edge_ids.clear();
    for (const int voxel : point_to_voxels[vi]) {
        for (const int id : voxel_to_primitives.at(voxel)) {
            if (is_edge_index(id)) {
                edge_ids.insert(to_edge_index(id));
            }
        }
    }
}

void SpatialHash::query_point_for_triangles(
    int vi, unordered_set<int>& tri_ids) const
{
    tri_ids.clear();
    for (const int voxel : point_to_voxels[vi]) {
        for (const int id : voxel_to_primitives.at(voxel)) {
            if (is_triangle_index(id)) {
                tri_ids.insert(to_triangle_index(id));
            }
        }
    }
}

// will only put edges with larger than eai index into edge_ids
void SpatialHash::query_edge_for_edges(
    int eai, unordered_set<int>& edge_ids) const
{
    edge_ids.clear();
    for (const int voxel : edge_to_voxels[eai]) {
        for (const int id : voxel_to_primitives.at(voxel)) {
            if (is_edge_index(id) && to_edge_index(id) > eai) {
                edge_ids.insert(to_edge_index(id));
            }
        }
    }
}

void SpatialHash::query_edge_for_triangles(
    int ei, unordered_set<int>& tri_ids) const
{
    tri_ids.clear();
    for (const int voxel : edge_to_voxels[ei]) {
        for (const auto& id : voxel_to_primitives.at(voxel)) {
            if (is_triangle_index(id)) {
                tri_ids.insert(to_triangle_index(id));
            }
        }
    }
}

// will only put triangles with larger than fai index into tri_ids
void SpatialHash::query_triangle_for_triangles(
    int fai, unordered_set<int>& tri_ids) const
{
    tri_ids.clear();
    for (const int voxel : face_to_voxels[fai]) {
        for (const auto& id : voxel_to_primitives.at(voxel)) {
            if (is_triangle_index(id) && to_triangle_index(id) > fai) {
                tri_ids.insert(to_triangle_index(id));
            }
        }
    }
}

// ============================================================================
// BroadPhase API

template <typename Candidate, bool swap_order, bool triangular>
void SpatialHash::detect_candidates(
    const std::vector<AABB>& boxesA,
    const std::vector<AABB>& boxesB,
    const std::function<void(int, unordered_set<int>&)>& query_A_for_Bs,
    const std::function<bool(int, int)>& can_collide,
    std::vector<Candidate>& candidates) const
{
    tbb::enumerable_thread_specific<std::vector<Candidate>> storage;

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), boxesA.size()),
        [&](const tbb::blocked_range<size_t>& range) {
            auto& local_candidates = storage.local();

            for (size_t i = range.begin(); i != range.end(); i++) {
                unordered_set<int> js;
                query_A_for_Bs(i, js);

                for (const int j : js) {
                    int ai = i, bi = j;
                    if constexpr (swap_order) {
                        std::swap(ai, bi);
                    }

                    if constexpr (triangular) {
                        if (ai >= bi) {
                            continue;
                        }
                    }

                    if (!can_collide(ai, bi)) {
                        continue;
                    }

                    if (boxesA[i].intersects(boxesB[j])) {
                        local_candidates.emplace_back(ai, bi);
                    }
                }
            }
        });

    merge_thread_local_vectors(storage, candidates);
}

template <typename Candidate>
void SpatialHash::detect_candidates(
    const std::vector<AABB>& boxesA,
    const std::function<void(int, unordered_set<int>&)>& query_A_for_As,
    const std::function<bool(int, int)>& can_collide,
    std::vector<Candidate>& candidates) const
{
    detect_candidates<Candidate, /*swap_order=*/false, /*triangular=*/true>(
        boxesA, boxesA, query_A_for_As, can_collide, candidates);
}

void SpatialHash::detect_vertex_vertex_candidates(
    std::vector<VertexVertexCandidate>& candidates) const
{
    if (vertex_boxes.size() == 0) {
        return;
    }

    detect_candidates(
        vertex_boxes,
        std::bind(&SpatialHash::query_point_for_points, this, _1, _2),
        can_vertices_collide, candidates);
}

void SpatialHash::detect_edge_vertex_candidates(
    std::vector<EdgeVertexCandidate>& candidates) const
{
    if (edge_boxes.size() == 0 || vertex_boxes.size() == 0) {
        return;
    }

    detect_candidates<EdgeVertexCandidate, /*swap_order=*/true>(
        vertex_boxes, edge_boxes,
        std::bind(&SpatialHash::query_point_for_edges, this, _1, _2),
        std::bind(&SpatialHash::can_edge_vertex_collide, this, _1, _2),
        candidates);
}

void SpatialHash::detect_edge_edge_candidates(
    std::vector<EdgeEdgeCandidate>& candidates) const
{
    if (edge_boxes.size() == 0) {
        return;
    }

    detect_candidates(
        edge_boxes, std::bind(&SpatialHash::query_edge_for_edges, this, _1, _2),
        std::bind(&SpatialHash::can_edges_collide, this, _1, _2), candidates);
}

void SpatialHash::detect_face_vertex_candidates(
    std::vector<FaceVertexCandidate>& candidates) const
{
    if (face_boxes.size() == 0 || vertex_boxes.size() == 0) {
        return;
    }

    // The ratio vertices:faces is 1:2, so we want to iterate over the vertices.
    detect_candidates<FaceVertexCandidate, /*swap_order=*/true>(
        vertex_boxes, face_boxes,
        std::bind(&SpatialHash::query_point_for_triangles, this, _1, _2),
        std::bind(&SpatialHash::can_face_vertex_collide, this, _1, _2),
        candidates);
}

void SpatialHash::detect_edge_face_candidates(
    std::vector<EdgeFaceCandidate>& candidates) const
{
    if (edge_boxes.size() == 0 || face_boxes.size() == 0) {
        return;
    }

    detect_candidates<EdgeFaceCandidate, /*swap_order=*/false>(
        edge_boxes, face_boxes,
        std::bind(&SpatialHash::query_edge_for_triangles, this, _1, _2),
        std::bind(&SpatialHash::can_edge_face_collide, this, _1, _2),
        candidates);
}

void SpatialHash::detect_face_face_candidates(
    std::vector<FaceFaceCandidate>& candidates) const
{
    if (face_boxes.size() == 0) {
        return;
    }

    detect_candidates(
        face_boxes,
        std::bind(&SpatialHash::query_triangle_for_triangles, this, _1, _2),
        std::bind(&SpatialHash::can_faces_collide, this, _1, _2), candidates);
}

// ============================================================================

int SpatialHash::locate_voxel_index(const VectorMax3d& p) const
{
    return voxel_axis_index_to_voxel_index(locate_voxel_axis_index(p));
}

ArrayMax3i SpatialHash::locate_voxel_axis_index(const VectorMax3d& p) const
{
    return ((p.array() - left_bottom_corner) * one_div_voxelSize)
        .floor()
        .template cast<int>();
}

void SpatialHash::locate_box_voxel_axis_index(
    ArrayMax3d min_corner, // input but will be modified
    ArrayMax3d max_corner, // input but will be modified
    ArrayMax3i& min_index, // output
    ArrayMax3i& max_index, // output
    const double inflation_radius) const
{
    AABB::conservative_inflation(min_corner, max_corner, inflation_radius);
    min_index = locate_voxel_axis_index(min_corner).max(ArrayMax3i::Zero(dim));
    max_index = locate_voxel_axis_index(max_corner).min(voxel_count - 1);
}

int SpatialHash::voxel_axis_index_to_voxel_index(
    const ArrayMax3i& voxel_axis_index) const
{
    return voxel_axis_index_to_voxel_index(
        voxel_axis_index[0], voxel_axis_index[1],
        voxel_axis_index.size() >= 3 ? voxel_axis_index[2] : 0);
}

int SpatialHash::voxel_axis_index_to_voxel_index(int ix, int iy, int iz) const
{
    assert(ix >= 0 && ix < voxel_count[0]);
    assert(iy >= 0 && iy < voxel_count[1]);
    assert(iz >= 0 && iz < (voxel_count.size() >= 3 ? voxel_count[2] : 1));
    return ix + iy * voxel_count[0] + iz * voxel_count_0x1;
}

} // namespace ipc
