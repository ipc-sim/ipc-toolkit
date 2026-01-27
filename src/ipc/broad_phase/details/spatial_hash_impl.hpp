#include <ipc/broad_phase/spatial_hash.hpp>
#include <ipc/utils/unordered_map_and_set.hpp>

namespace ipc {

struct SpatialHash::Impl {
    // --- Data ---------------------------------------------------------------

    /// @brief Map from voxel index to the primitive indices it contains.
    unordered_map<int, std::vector<int>> voxel_to_primitives;

    /// @brief Map from point index to the voxel indices it occupies.
    std::vector<std::vector<int>> point_to_voxels;

    /// @brief Map from edge index to the voxel indices it occupies.
    std::vector<std::vector<int>> edge_to_voxels;

    /// @brief Map from face index to the voxel indices it occupies.
    std::vector<std::vector<int>> face_to_voxels;

    /// @brief The index of the first edge in voxel_occupancies
    int edge_start_ind = -1;

    /// @brief The index of the first triangle in voxel_occupancies
    int tri_start_ind = -1;

    // --- Methods ------------------------------------------------------------

    void fill_voxel_to_primitives()
    {
        edge_start_ind = point_to_voxels.size();
        tri_start_ind = edge_start_ind + edge_to_voxels.size();

        for (size_t i = 0; i < point_to_voxels.size(); i++) {
            for (const auto& voxel : point_to_voxels[i]) {
                voxel_to_primitives[voxel].emplace_back(i);
            }
        }

        for (size_t i = 0; i < edge_to_voxels.size(); i++) {
            for (const auto& voxel : edge_to_voxels[i]) {
                voxel_to_primitives[voxel].emplace_back(i + edge_start_ind);
            }
        }

        for (size_t i = 0; i < face_to_voxels.size(); i++) {
            for (const auto& voxel : face_to_voxels[i]) {
                voxel_to_primitives[voxel].emplace_back(i + tri_start_ind);
            }
        }
    }

    void clear()
    {
        voxel_to_primitives.clear();
        point_to_voxels.clear();
        edge_to_voxels.clear();
        face_to_voxels.clear();
    }

    /// @brief Check if primitive index refers to a vertex.
    bool is_vertex_index(int idx) const { return idx < edge_start_ind; }

    /// @brief Check if primitive index refers to an edge.
    bool is_edge_index(int idx) const
    {
        return idx >= edge_start_ind && idx < tri_start_ind;
    }

    /// @brief Check if primitive index refers to a triangle.
    bool is_triangle_index(int idx) const { return idx >= tri_start_ind; }

    /// @brief Convert a primitive index to an edge index.
    int to_edge_index(int idx) const
    {
        assert(is_edge_index(idx));
        return idx - edge_start_ind;
    }

    /// @brief Convert a primitive index to a triangle index.
    int to_triangle_index(int idx) const
    {
        assert(is_triangle_index(idx));
        return idx - tri_start_ind;
    }

    void query_point_for_points(int vi, unordered_set<int>& vert_ids) const
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

    void query_point_for_edges(int vi, unordered_set<int>& edge_ids) const
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

    void query_point_for_triangles(int vi, unordered_set<int>& tri_ids) const
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

    // will only put edges with larger than ei index into edge_ids
    void query_edge_for_edges(int eai, unordered_set<int>& edge_ids) const
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

    void query_edge_for_triangles(int ei, unordered_set<int>& tri_ids) const
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

    // will only put triangles with larger than ti index into tri_ids
    void query_triangle_for_triangles(int ti, unordered_set<int>& tri_ids) const
    {
        tri_ids.clear();
        for (const int voxel : face_to_voxels[ti]) {
            for (const auto& id : voxel_to_primitives.at(voxel)) {
                if (is_triangle_index(id) && to_triangle_index(id) > ti) {
                    tri_ids.insert(to_triangle_index(id));
                }
            }
        }
    }
};

} // namespace ipc
