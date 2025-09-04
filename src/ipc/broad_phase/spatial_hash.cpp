// Modified version of SpatialHash.hpp from IPC codebase.
// Originally created by Minchen Li.
#include "spatial_hash.hpp"

#include "details/spatial_hash_impl.hpp"

#include <ipc/config.hpp>
#include <ipc/broad_phase/voxel_size_heuristic.hpp>
#include <ipc/ccd/aabb.hpp>
#include <ipc/utils/merge_thread_local.hpp>

#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>

using namespace std::placeholders;

namespace ipc {

SpatialHash::SpatialHash() : impl(std::make_unique<SpatialHash::Impl>()) { }

SpatialHash::~SpatialHash() = default;

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
        Eigen::ConstRef<Eigen::Array3i> min_voxel,
        Eigen::ConstRef<Eigen::Array3i> max_voxel,
        Eigen::ConstRef<ArrayMax3i> voxel_count,
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
} // namespace

void SpatialHash::build(
    Eigen::ConstRef<Eigen::MatrixXd> vertices,
    Eigen::ConstRef<Eigen::MatrixXi> edges,
    Eigen::ConstRef<Eigen::MatrixXi> faces,
    double inflation_radius,
    double voxel_size)
{
    clear();
    build_vertex_boxes(vertices, vertex_boxes, inflation_radius);
    build(edges, faces, voxel_size);
}

void SpatialHash::build(
    Eigen::ConstRef<Eigen::MatrixXd> vertices_t0,
    Eigen::ConstRef<Eigen::MatrixXd> vertices_t1,
    Eigen::ConstRef<Eigen::MatrixXi> edges,
    Eigen::ConstRef<Eigen::MatrixXi> faces,
    double inflation_radius,
    double voxel_size)
{
    clear();
    build_vertex_boxes(
        vertices_t0, vertices_t1, vertex_boxes, inflation_radius);
    build(edges, faces, voxel_size);
}

void SpatialHash::build(
    const std::vector<AABB>& _vertex_boxes,
    Eigen::ConstRef<Eigen::MatrixXi> edges,
    Eigen::ConstRef<Eigen::MatrixXi> faces,
    double voxel_size)
{
    // WARNING: Clear will reset vertex_boxes if this assert is triggered
    assert(&(this->vertex_boxes) != &_vertex_boxes);
    clear();
    this->vertex_boxes = _vertex_boxes;
    build(edges, faces, voxel_size);
}

void SpatialHash::build(
    Eigen::ConstRef<Eigen::MatrixXi> edges,
    Eigen::ConstRef<Eigen::MatrixXi> faces,
    double voxel_size)
{
    const size_t num_vertices = vertex_boxes.size();
    assert(num_vertices > 0);
    dim = vertex_boxes[0].min.size();

    BroadPhase::build(edges, faces);

    if (voxel_size <= 0) {
        voxel_size = suggest_good_voxel_size(
            edges.rows() > 0 ? edge_boxes : vertex_boxes);
    }

    left_bottom_corner = vertex_boxes[0].min;
    right_top_corner = vertex_boxes[0].max;
    for (const auto& box : vertex_boxes) {
        left_bottom_corner = left_bottom_corner.min(box.min);
        right_top_corner = right_top_corner.max(box.max);
    }

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
        vertex_min_voxel_axis_index[vi].head(dim) =
            locate_voxel_axis_index(vertex_boxes[vi].min);
        vertex_max_voxel_axis_index[vi].head(dim) =
            locate_voxel_axis_index(vertex_boxes[vi].max);
    });

    // ------------------------------------------------------------------------

    impl->point_to_voxels.resize(num_vertices);
    tbb_parallel_block_range_for(0, num_vertices, [&](size_t vi) {
        fill_primitive_to_voxels(
            vertex_min_voxel_axis_index[vi], vertex_max_voxel_axis_index[vi],
            voxel_count, voxel_count_0x1, impl->point_to_voxels[vi]);
    });

    impl->edge_to_voxels.resize(edges.rows());
    tbb_parallel_block_range_for(0, edges.rows(), [&](size_t ei) {
        fill_primitive_to_voxels(
            vertex_min_voxel_axis_index[edges(ei, 0)].min(
                vertex_min_voxel_axis_index[edges(ei, 1)]),
            vertex_max_voxel_axis_index[edges(ei, 0)].max(
                vertex_max_voxel_axis_index[edges(ei, 1)]),
            voxel_count, voxel_count_0x1, impl->edge_to_voxels[ei]);
    });

    impl->face_to_voxels.resize(faces.rows());
    tbb_parallel_block_range_for(0, faces.rows(), [&](size_t fi) {
        fill_primitive_to_voxels(
            vertex_min_voxel_axis_index[faces(fi, 0)]
                .min(vertex_min_voxel_axis_index[faces(fi, 1)])
                .min(vertex_min_voxel_axis_index[faces(fi, 2)]),
            vertex_max_voxel_axis_index[faces(fi, 0)]
                .max(vertex_max_voxel_axis_index[faces(fi, 1)])
                .max(vertex_max_voxel_axis_index[faces(fi, 2)]),
            voxel_count, voxel_count_0x1, impl->face_to_voxels[fi]);
    });

    // ------------------------------------------------------------------------

    impl->fill_voxel_to_primitives();
}

void SpatialHash::clear()
{
    BroadPhase::clear();
    impl->clear();
}

// ============================================================================
// BroadPhase API

namespace {

    template <typename Candidate, bool swap_order, bool triangular = false>
    void detect_candidates(
        const std::vector<AABB>& boxesA,
        const std::vector<AABB>& boxesB,
        const std::function<void(int, unordered_set<int>&)>& query_A_for_Bs,
        const std::function<bool(size_t, size_t)>& can_collide,
        std::vector<Candidate>& candidates)
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
    void detect_candidates(
        const std::vector<AABB>& boxesA,
        const std::function<void(int, unordered_set<int>&)>& query_A_for_As,
        const std::function<bool(size_t, size_t)>& can_collide,
        std::vector<Candidate>& candidates)
    {
        detect_candidates<Candidate, /*swap_order=*/false, /*triangular=*/true>(
            boxesA, boxesA, query_A_for_As, can_collide, candidates);
    }

} // namespace

void SpatialHash::detect_vertex_vertex_candidates(
    std::vector<VertexVertexCandidate>& candidates) const
{
    if (vertex_boxes.size() == 0) {
        return;
    }

    detect_candidates(
        vertex_boxes,
        [&](int i, unordered_set<int>& ids) {
            impl->query_point_for_points(i, ids);
        },
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
        [&](int i, unordered_set<int>& ids) {
            impl->query_point_for_edges(i, ids);
        },
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
        edge_boxes,
        [&](int i, unordered_set<int>& ids) {
            impl->query_edge_for_edges(i, ids);
        },
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
        [&](int i, unordered_set<int>& ids) {
            impl->query_point_for_triangles(i, ids);
        },
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
        [&](int i, unordered_set<int>& ids) {
            impl->query_edge_for_triangles(i, ids);
        },
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
        [&](int i, unordered_set<int>& ids) {
            impl->query_triangle_for_triangles(i, ids);
        },
        std::bind(&SpatialHash::can_faces_collide, this, _1, _2), candidates);
}

// ============================================================================

int SpatialHash::locate_voxel_index(Eigen::ConstRef<VectorMax3d> p) const
{
    return voxel_axis_index_to_voxel_index(locate_voxel_axis_index(p));
}

ArrayMax3i
SpatialHash::locate_voxel_axis_index(Eigen::ConstRef<VectorMax3d> p) const
{
    return ((p.array() - left_bottom_corner) * one_div_voxelSize)
        .floor()
        .template cast<int>();
}

void SpatialHash::locate_box_voxel_axis_index(
    ArrayMax3d min_corner,            // input but will be modified
    ArrayMax3d max_corner,            // input but will be modified
    Eigen::Ref<ArrayMax3i> min_index, // output
    Eigen::Ref<ArrayMax3i> max_index, // output
    const double inflation_radius) const
{
    AABB::conservative_inflation(min_corner, max_corner, inflation_radius);
    min_index = locate_voxel_axis_index(min_corner).max(ArrayMax3i::Zero(dim));
    max_index = locate_voxel_axis_index(max_corner).min(voxel_count - 1);
}

int SpatialHash::voxel_axis_index_to_voxel_index(
    Eigen::ConstRef<ArrayMax3i> voxel_axis_index) const
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
