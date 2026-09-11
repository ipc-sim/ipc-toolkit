#pragma once

#include <ipc/collision_mesh.hpp>
#include <ipc/broad_phase/lbvh.hpp>

#include <vector>

namespace ipc {

/// @brief Broad-phase index over a mesh's vertices/edges/faces, built once
/// per vertex configuration and queried many times for arbitrary points in
/// space.
///
/// Wraps ipc::LBVH, reusing its Node/Nodes data structures and build() as
/// they are; adds the point-vs-tree query LBVH's public API does not expose
/// (LBVH's own traversal is a file-local template, not part of its public
/// interface).
///
/// The tree is built with zero box inflation (mirroring the pattern used by
/// the potential-sampling tools this class is ported from): the query
/// radius is applied to the query point's own box instead. This decouples
/// the built tree from any particular dhat, so one build serves queries at
/// any radius.
class ArbitraryPointBVH {
public:
    /// @brief Rebuild the tree. O(n log n). Call once per vertex
    /// configuration (i.e. whenever V changes).
    void update(Eigen::ConstRef<Eigen::MatrixXd> V, const CollisionMesh& mesh);

    /// @brief Find real mesh primitives whose (uninflated) AABB is within
    /// `radius` of q, i.e. intersects a box of half-width `radius` centered
    /// at q. Broad-phase only: callers still need to check exact distances.
    /// @param q Query point (2D or 3D, matching the mesh).
    /// @param radius Half-width of the query box around q.
    /// @param[out] vertex_ids Real vertex ids found nearby.
    /// @param[out] edge_ids Real edge ids found nearby.
    /// @param[out] face_ids Real face ids found nearby. Always empty in 2D,
    ///   where mesh.faces() is empty and no face tree is built.
    void query_point(
        Eigen::ConstRef<RowVectorMax3d> q,
        double radius,
        std::vector<index_t>& vertex_ids,
        std::vector<index_t>& edge_ids,
        std::vector<index_t>& face_ids) const;

private:
    LBVH bvh;
};

} // namespace ipc
