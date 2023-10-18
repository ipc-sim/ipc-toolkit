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

namespace ipc {

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
    ArrayMax3d range = right_top_corner - left_bottom_corner;
    voxel_count = (range * one_div_voxelSize).ceil().template cast<int>();
    if (voxel_count.minCoeff() <= 0) {
        // cast overflow due to huge search direction
        one_div_voxelSize = 1.0 / (range.maxCoeff() * 1.01);
        voxel_count.setOnes();
    }
    voxel_count_0x1 = voxel_count[0] * voxel_count[1];

    edge_start_ind = num_vertices;
    tri_start_ind = edge_start_ind + edges.rows();

    // precompute vVAI
    std::vector<Eigen::Array3i> vertexMinVAI(
        num_vertices, Eigen::Array3i::Zero());
    std::vector<Eigen::Array3i> vertexMaxVAI(
        num_vertices, Eigen::Array3i::Zero());
    tbb::parallel_for(size_t(0), num_vertices, [&](size_t vi) {
        ArrayMax3d v_min = vertices_t0.row(vi).cwiseMin(vertices_t1.row(vi));
        ArrayMax3d v_max = vertices_t0.row(vi).cwiseMax(vertices_t1.row(vi));
        AABB::conservative_inflation(v_min, v_max, inflation_radius);

        ArrayMax3i v_vai_min, v_vai_max;
        locate_voxel_axis_index(v_min, v_vai_min);
        locate_voxel_axis_index(v_max, v_vai_max);

        vertexMinVAI[vi].head(dim) = v_vai_min;
        vertexMaxVAI[vi].head(dim) = v_vai_max;
    });

    point_and_edge_occupancy.resize(tri_start_ind);

    tbb::parallel_for(size_t(0), num_vertices, [&](size_t vi) {
        const Eigen::Array3i &mins = vertexMinVAI[vi], &maxs = vertexMaxVAI[vi];
        assert((mins <= maxs).all());
        point_and_edge_occupancy[vi].reserve((maxs - mins + 1).prod());
        for (int iz = mins[2]; iz <= maxs[2]; iz++) {
            int z_offset = iz * voxel_count_0x1;
            for (int iy = mins[1]; iy <= maxs[1]; iy++) {
                int yz_offset = iy * voxel_count[0] + z_offset;
                for (int ix = mins[0]; ix <= maxs[0]; ix++) {
                    point_and_edge_occupancy[vi].emplace_back(ix + yz_offset);
                }
            }
        }
    });

    tbb::parallel_for(size_t(0), size_t(edges.rows()), [&](size_t ei) {
        int ei_ind = ei + edge_start_ind;

        Eigen::Array3i mins =
            vertexMinVAI[edges(ei, 0)].min(vertexMinVAI[edges(ei, 1)]);
        Eigen::Array3i maxs =
            vertexMaxVAI[edges(ei, 0)].max(vertexMaxVAI[edges(ei, 1)]);

        point_and_edge_occupancy[ei_ind].reserve((maxs - mins + 1).prod());
        for (int iz = mins[2]; iz <= maxs[2]; iz++) {
            int z_offset = iz * voxel_count_0x1;
            for (int iy = mins[1]; iy <= maxs[1]; iy++) {
                int yz_offset = iy * voxel_count[0] + z_offset;
                for (int ix = mins[0]; ix <= maxs[0]; ix++) {
                    point_and_edge_occupancy[ei_ind].emplace_back(
                        ix + yz_offset);
                }
            }
        }
    });

    std::vector<std::vector<int>> voxelLoc_f(faces.rows());
    tbb::parallel_for(size_t(0), size_t(faces.rows()), [&](size_t fi) {
        Eigen::Array3i mins = vertexMinVAI[faces(fi, 0)]
                                  .min(vertexMinVAI[faces(fi, 1)])
                                  .min(vertexMinVAI[faces(fi, 2)]);
        Eigen::Array3i maxs = vertexMaxVAI[faces(fi, 0)]
                                  .max(vertexMaxVAI[faces(fi, 1)])
                                  .max(vertexMaxVAI[faces(fi, 2)]);

        for (int iz = mins[2]; iz <= maxs[2]; iz++) {
            int z_offset = iz * voxel_count_0x1;
            for (int iy = mins[1]; iy <= maxs[1]; iy++) {
                int yz_offset = iy * voxel_count[0] + z_offset;
                for (int ix = mins[0]; ix <= maxs[0]; ix++) {
                    voxelLoc_f[fi].emplace_back(ix + yz_offset);
                }
            }
        }
    });

    for (int i = 0; i < point_and_edge_occupancy.size(); i++) {
        for (const auto& voxelI : point_and_edge_occupancy[i]) {
            voxel[voxelI].emplace_back(i);
        }
    }
    for (int fi = 0; fi < voxelLoc_f.size(); fi++) {
        for (const auto& voxelI : voxelLoc_f[fi]) {
            voxel[voxelI].emplace_back(fi + tri_start_ind);
        }
    }
}

void SpatialHash::query_point_for_points(
    int vi, unordered_set<int>& vert_inds) const
{
    vert_inds.clear();
    for (const auto& voxel_ind : point_and_edge_occupancy[vi]) {
        const auto& voxelI = voxel.at(voxel_ind);
        for (const auto& indI : voxelI) {
            if (is_vertex_index(indI) && indI > vi) {
                vert_inds.insert(indI);
            }
        }
    }
}

void SpatialHash::query_point_for_edges(
    int vi, unordered_set<int>& edge_inds) const
{
    edge_inds.clear();
    for (const auto& voxel_ind : point_and_edge_occupancy[vi]) {
        const auto& voxelI = voxel.at(voxel_ind);
        for (const auto& indI : voxelI) {
            if (is_edge_index(indI)) {
                edge_inds.insert(to_edge_index(indI));
            }
        }
    }
}

void SpatialHash::query_point_for_triangles(
    int vi, unordered_set<int>& tri_inds) const
{
    tri_inds.clear();
    for (const auto& voxel_ind : point_and_edge_occupancy[vi]) {
        const auto& voxelI = voxel.at(voxel_ind);
        for (const auto& indI : voxelI) {
            if (is_triangle_index(indI)) {
                tri_inds.insert(to_triangle_index(indI));
            }
        }
    }
}

// will only put edges with larger than eai index into edge_inds
void SpatialHash::query_edge_for_edges(
    int eai, unordered_set<int>& edge_inds) const
{
    edge_inds.clear();
    for (const auto& voxel_ind :
         point_and_edge_occupancy[eai + edge_start_ind]) {
        const auto& voxelI = voxel.at(voxel_ind);
        for (const auto& indI : voxelI) {
            if (is_edge_index(indI) && to_edge_index(indI) > eai) {
                edge_inds.insert(to_edge_index(indI));
            }
        }
    }
}

void SpatialHash::query_edge_for_triangles(
    int ei, unordered_set<int>& tri_inds) const
{
    tri_inds.clear();
    for (const auto& voxel_ind :
         point_and_edge_occupancy[ei + edge_start_ind]) {
        const auto& voxelI = voxel.at(voxel_ind);
        for (const auto& indI : voxelI) {
            if (is_triangle_index(indI)) {
                tri_inds.insert(to_triangle_index(indI));
            }
        }
    }
}

// ============================================================================
// BroadPhase API

void SpatialHash::detect_vertex_vertex_candidates(
    std::vector<VertexVertexCandidate>& candidates) const
{
    tbb::enumerable_thread_specific<std::vector<VertexVertexCandidate>>
        storages;

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), vertex_boxes.size()),
        [&](const tbb::blocked_range<size_t>& range) {
            auto& local_candidates = storages.local();

            for (long vi = range.begin(); vi != range.end(); vi++) {
                unordered_set<int> vertex_inds;
                query_point_for_points(vi, vertex_inds);

                for (const auto& vj : vertex_inds) {
                    if (vi >= vj || !can_vertices_collide(vi, vj)) {
                        continue;
                    }

                    if (vertex_boxes[vi].intersects(vertex_boxes[vi])) {
                        local_candidates.emplace_back(vi, vj);
                    }
                }
            }
        });

    merge_thread_local_vectors(storages, candidates);
}

void SpatialHash::detect_edge_vertex_candidates(
    std::vector<EdgeVertexCandidate>& candidates) const
{
    tbb::enumerable_thread_specific<std::vector<EdgeVertexCandidate>> storages;

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), vertex_boxes.size()),
        [&](const tbb::blocked_range<size_t>& range) {
            auto& local_candidates = storages.local();

            for (long vi = range.begin(); vi != range.end(); vi++) {
                const AABB& vertex_box = vertex_boxes[vi];

                unordered_set<int> edge_inds;
                query_point_for_edges(vi, edge_inds);

                for (const auto& ei : edge_inds) {
                    if (!can_edge_vertex_collide(ei, vi)) {
                        continue;
                    }

                    const AABB& edge_box = edge_boxes[ei];
                    if (vertex_box.intersects(edge_box)) {
                        local_candidates.emplace_back(ei, vi);
                    }
                }
            }
        });

    merge_thread_local_vectors(storages, candidates);
}

void SpatialHash::detect_edge_edge_candidates(
    std::vector<EdgeEdgeCandidate>& candidates) const
{
    tbb::enumerable_thread_specific<std::vector<EdgeEdgeCandidate>> storages;

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), edge_boxes.size()),
        [&](const tbb::blocked_range<size_t>& range) {
            auto& local_candidates = storages.local();

            for (long eai = range.begin(); eai != range.end(); eai++) {
                const AABB& edge_a_box = edge_boxes[eai];

                unordered_set<int> edge_inds;
                query_edge_for_edges(eai, edge_inds);

                for (const auto& ebi : edge_inds) {
                    if (eai >= ebi || !can_edges_collide(eai, ebi)) {
                        continue;
                    }

                    const AABB& edge_b_box = edge_boxes[ebi];
                    if (edge_a_box.intersects(edge_b_box)) {
                        local_candidates.emplace_back(eai, ebi);
                    }
                }
            }
        });

    merge_thread_local_vectors(storages, candidates);
}

void SpatialHash::detect_face_vertex_candidates(
    std::vector<FaceVertexCandidate>& candidates) const
{
    tbb::enumerable_thread_specific<std::vector<FaceVertexCandidate>> storages;

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), vertex_boxes.size()),
        [&](const tbb::blocked_range<size_t>& range) {
            auto& local_candidates = storages.local();

            for (long vi = range.begin(); vi != range.end(); vi++) {
                const AABB& vertex_box = vertex_boxes[vi];

                unordered_set<int> tri_inds;
                query_point_for_triangles(vi, tri_inds);

                for (const auto& fi : tri_inds) {
                    if (!can_face_vertex_collide(fi, vi)) {
                        continue;
                    }

                    const AABB& face_box = face_boxes[fi];
                    if (vertex_box.intersects(face_box)) {
                        local_candidates.emplace_back(fi, vi);
                    }
                }
            }
        });

    merge_thread_local_vectors(storages, candidates);
}

void SpatialHash::detect_edge_face_candidates(
    std::vector<EdgeFaceCandidate>& candidates) const
{
    tbb::enumerable_thread_specific<std::vector<EdgeFaceCandidate>> storages;

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), edge_boxes.size()),
        [&](const tbb::blocked_range<size_t>& range) {
            auto& local_candidates = storages.local();

            for (long ei = range.begin(); ei != range.end(); ei++) {
                const AABB& edge_box = edge_boxes[ei];

                unordered_set<int> tri_inds;
                query_edge_for_triangles(ei, tri_inds);

                for (const auto& fi : tri_inds) {
                    if (!can_edge_face_collide(ei, fi)) {
                        continue;
                    }

                    const AABB& face_box = face_boxes[fi];
                    if (edge_box.intersects(face_box)) {
                        local_candidates.emplace_back(ei, fi);
                    }
                }
            }
        });

    merge_thread_local_vectors(storages, candidates);
}

// ============================================================================

int SpatialHash::locate_voxel_index(const VectorMax3d& p) const
{
    ArrayMax3i voxel_axis_index;
    locate_voxel_axis_index(p, voxel_axis_index);
    return voxelAxisIndex2VoxelIndex(voxel_axis_index);
}

void SpatialHash::SpatialHash::locate_voxel_axis_index(
    const VectorMax3d& p, ArrayMax3i& voxel_axis_index) const
{
    voxel_axis_index = ((p.array() - left_bottom_corner) * one_div_voxelSize)
                           .floor()
                           .template cast<int>();
}

void SpatialHash::locate_box_voxel_axis_index(
    ArrayMax3d min_corner,
    ArrayMax3d max_corner,
    ArrayMax3i& min_index,
    ArrayMax3i& max_index,
    const double inflation_radius) const
{
    AABB::conservative_inflation(min_corner, max_corner, inflation_radius);
    locate_voxel_axis_index(min_corner, min_index);
    locate_voxel_axis_index(max_corner, max_index);
    min_index = min_index.max(ArrayMax3i::Zero(dim));
    max_index = max_index.min(voxel_count - 1);
}

int SpatialHash::voxelAxisIndex2VoxelIndex(
    const ArrayMax3i& voxel_axis_index) const
{
    return voxelAxisIndex2VoxelIndex(
        voxel_axis_index[0], voxel_axis_index[1],
        voxel_axis_index.size() >= 3 ? voxel_axis_index[2] : 0);
}

int SpatialHash::voxelAxisIndex2VoxelIndex(int ix, int iy, int iz) const
{
    assert(ix >= 0 && ix < voxel_count[0]);
    assert(iy >= 0 && iy < voxel_count[1]);
    assert(iz >= 0 && iz < (voxel_count.size() >= 3 ? voxel_count[2] : 1));
    return ix + iy * voxel_count[0] + iz * voxel_count_0x1;
}

} // namespace ipc
