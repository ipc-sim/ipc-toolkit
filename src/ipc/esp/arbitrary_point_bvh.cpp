#include "arbitrary_point_bvh.hpp"

#include <ipc/broad_phase/details/lbvh_traverse.hpp>

#include <vector>

namespace ipc {

namespace {

    /// @brief Build a query "node" whose AABB is a cube of half-width
    /// `radius` centered at p, for use with LBVH::Node::intersects().
    /// Only aabb_min/aabb_max are meaningful; the union fields are unused
    /// by intersects() but are set to keep the node well-defined.
    ///
    /// A 2D point is padded to (x, y, 0), which is what AABB does to 2D
    /// primitives (it stores Array3d and zero-fills the unused component),
    /// so the resulting z range [-radius, radius] straddles the tree's.
    LBVH::Node
    point_query_node(Eigen::ConstRef<RowVectorMax3d> p, double radius)
    {
        Eigen::Array3d c = Eigen::Array3d::Zero();
        c.head(p.size()) = p.transpose().array();

        LBVH::Node n;
        n.aabb_min = (c - radius).cast<float>();
        n.aabb_max = (c + radius).cast<float>();
        n.primitive_id = -1;
        n.is_inner_marker = 0;
        return n;
    }

    /// @brief Collect the primitive ids of every leaf of one LBVH tree whose
    /// AABB intersects the query box. We reuse LBVH's own descent
    /// (details::traverse_lbvh), whose fixed-size stack is sized for the
    /// deepest tree a Morton-code build can produce and throws instead of
    /// overflowing.
    void query_point_vs_bvh(
        const LBVH::Nodes& bvh,
        const LBVH::Node& query,
        std::vector<index_t>& hits)
    {
        if (bvh.empty()) {
            return;
        }

        details::traverse_lbvh</*triangular=*/false>(
            /*query_leaf_idx=*/-1, bvh.data(), static_cast<int>(bvh.size()),
            /*target_rightmost=*/nullptr,
            [&](const LBVH::Node& node) { return node.intersects(query); },
            [&](const LBVH::Node& leaf, int /*leaf_idx*/, bool /*mask*/) {
                hits.push_back(leaf.primitive_id);
            });
    }

} // namespace

void ArbitraryPointBVH::update(
    Eigen::ConstRef<Eigen::MatrixXd> V, const CollisionMesh& mesh)
{
    // Zero inflation: the query radius is applied to the query point's own
    // box at query time instead, so a single build serves queries at any
    // radius.
    m_bvh.build(V, mesh.edges(), mesh.faces(), /*inflation_radius=*/0.0);
}

void ArbitraryPointBVH::query_point(
    Eigen::ConstRef<RowVectorMax3d> q,
    double radius,
    std::vector<index_t>& vertex_ids,
    std::vector<index_t>& edge_ids,
    std::vector<index_t>& face_ids) const
{
    vertex_ids.clear();
    edge_ids.clear();
    face_ids.clear();

    const LBVH::Node qnode = point_query_node(q, radius);
    query_point_vs_bvh(m_bvh.vertex_nodes(), qnode, vertex_ids);
    query_point_vs_bvh(m_bvh.edge_nodes(), qnode, edge_ids);
    query_point_vs_bvh(m_bvh.face_nodes(), qnode, face_ids);
}

} // namespace ipc
