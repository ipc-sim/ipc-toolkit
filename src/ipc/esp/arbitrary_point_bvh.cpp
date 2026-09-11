#include "arbitrary_point_bvh.hpp"

#include <cassert>
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

    /// @brief Stack-based traversal of a single LBVH tree against one query
    /// box, collecting the primitive ids of every leaf whose AABB
    /// intersects the query. LBVH's own traversal (lbvh.cpp) is a
    /// file-local template used for tree-vs-tree queries; this is the
    /// tree-vs-single-external-box equivalent, which LBVH does not expose.
    void query_point_vs_bvh(
        const LBVH::Nodes& bvh,
        const LBVH::Node& query,
        std::vector<index_t>& hits)
    {
        if (bvh.empty()) {
            return;
        }

        // A single-primitive tree is one node, and that node is a leaf, so
        // the descend-from-an-inner-root loop below would read the union's
        // left/right through a leaf. Easy to hit in 2D (a one-edge mesh).
        if (bvh.size() == 1) {
            if (bvh[0].intersects(query)) {
                hits.push_back(bvh[0].primitive_id);
            }
            return;
        }

        constexpr int MAX_STACK = 64;
        int stack[MAX_STACK];
        int sp = 0;
        stack[sp++] = LBVH::Node::INVALID_POINTER;

        int node_idx = 0; // root is always index 0
        do {
            const LBVH::Node& node = bvh[node_idx];
            assert(node.is_inner());

            const LBVH::Node& cl = bvh[node.left];
            const LBVH::Node& cr = bvh[node.right];
            const bool hl = cl.intersects(query);
            const bool hr = cr.intersects(query);

            if (hl && cl.is_leaf()) {
                hits.push_back(cl.primitive_id);
            }
            if (hr && cr.is_leaf()) {
                hits.push_back(cr.primitive_id);
            }

            const bool tl = hl && !cl.is_leaf();
            const bool tr = hr && !cr.is_leaf();

            if (!tl && !tr) {
                node_idx = stack[--sp];
            } else {
                node_idx = tl ? node.left : node.right;
                if (tl && tr) {
                    assert(sp < MAX_STACK);
                    stack[sp++] = node.right;
                }
            }
        } while (node_idx != LBVH::Node::INVALID_POINTER);
    }

} // namespace

void ArbitraryPointBVH::update(
    Eigen::ConstRef<Eigen::MatrixXd> V, const CollisionMesh& mesh)
{
    // Zero inflation: the query radius is applied to the query point's own
    // box at query time instead, so a single build serves queries at any
    // radius.
    bvh.build(V, mesh.edges(), mesh.faces(), /*inflation_radius=*/0.0);
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
    query_point_vs_bvh(bvh.vertex_nodes(), qnode, vertex_ids);
    query_point_vs_bvh(bvh.edge_nodes(), qnode, edge_ids);
    query_point_vs_bvh(bvh.face_nodes(), qnode, face_ids);
}

} // namespace ipc
