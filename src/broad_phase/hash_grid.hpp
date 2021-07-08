#pragma once

#include <vector>

#include <Eigen/Core>

#include <ipc/broad_phase/collision_candidate.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

/// @brief Axis aligned bounding-box of some type
class AABB {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    AABB() {}

    AABB(const ArrayMax3d& min, const ArrayMax3d& max);

    AABB(const AABB& aabb1, const AABB& aabb2)
        : AABB(aabb1.min.min(aabb2.min), aabb1.max.max(aabb2.max))
    {
    }

    AABB(const AABB& aabb1, const AABB& aabb2, const AABB& aabb3)
        : AABB(
              aabb1.min.min(aabb2.min).min(aabb3.min),
              aabb1.max.max(aabb2.max).max(aabb3.max))
    {
    }

    static bool are_overlapping(const AABB& a, const AABB& b);

    inline const ArrayMax3d& getMin() const { return min; }
    inline const ArrayMax3d& getMax() const { return max; }
    inline const ArrayMax3d& getHalfExtent() const { return half_extent; }
    inline const ArrayMax3d& getCenter() const { return center; }

protected:
    ArrayMax3d min;
    ArrayMax3d max;
    ArrayMax3d half_extent;
    ArrayMax3d center;
};

/// @brief An entry into the hash grid as a (key, value) pair.
class HashItem {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    int key;   /// @brief The key of the item.
    int id;    /// @brief The value of the item.
    AABB aabb; /// @brief The axis-aligned bounding box of the element

    /// @brief Construct a hash item as a (key, value) pair.
    HashItem(int key, int id, const AABB& aabb)
        : key(key)
        , id(id)
        , aabb(aabb)
    {
    }

    /// @brief Compare HashItems by their keys for sorting.
    bool operator<(const HashItem& other) const
    {
        if (key == other.key) {
            return id < other.id;
        }
        return key < other.key;
    }
};

class HashGrid {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    double cellSize() const { return m_cellSize; }
    const ArrayMax3i& gridSize() const { return m_gridSize; }
    const ArrayMax3d& domainMin() const { return m_domainMin; }
    const ArrayMax3d& domainMax() const { return m_domainMax; }

    void resizeFromBox(
        const ArrayMax3d& min, const ArrayMax3d& max, double cellSize);

    void resize(
        const Eigen::MatrixXd& vertices_t0,
        const Eigen::MatrixXd& vertices_t1,
        const Eigen::MatrixXi& edges,
        const double inflation_radius = 0.0);

    void resize(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const double inflation_radius = 0.0)
    {
        resize(vertices, vertices, edges, inflation_radius);
    }

    /// @brief Add a vertex as a AABB containing the time swept edge.
    void addVertex(
        const VectorMax3d& vertex_t0,
        const VectorMax3d& vertex_t1,
        const long index,
        const double inflation_radius = 0.0);

    /// @brief Add a vertex.
    void addVertex(
        const VectorMax3d& vertex,
        const long index,
        const double inflation_radius = 0.0)
    {
        addVertex(vertex, vertex, index, inflation_radius);
    }

    /// @brief Add all vertices as AABBs containing the time swept edge.
    void addVertices(
        const Eigen::MatrixXd& vertices_t0,
        const Eigen::MatrixXd& vertices_t1,
        const double inflation_radius = 0.0);

    /// @brief Add all vertices.
    void addVertices(
        const Eigen::MatrixXd& vertices, const double inflation_radius = 0.0)
    {
        addVertices(vertices, vertices, inflation_radius);
    }

    /// @brief Add all vertices as AABBs containing the time swept edge.
    void addVerticesFromEdges(
        const Eigen::MatrixXd& vertices_t0,
        const Eigen::MatrixXd& vertices_t1,
        const Eigen::MatrixXi& edges,
        const double inflation_radius = 0.0);

    /// @brief Add all vertices as AABBs containing the time swept edge.
    void addVerticesFromEdges(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const double inflation_radius = 0.0)
    {
        addVerticesFromEdges(vertices, vertices, edges, inflation_radius);
    }

    /// @brief Add an edge as a AABB containing the time swept quad.
    void addEdge(
        const VectorMax3d& edge_vertex0_t0,
        const VectorMax3d& edge_vertex1_t0,
        const VectorMax3d& edge_vertex0_t1,
        const VectorMax3d& edge_vertex1_t1,
        const long index,
        const double inflation_radius = 0.0);

    /// @brief Add an edge as a AABB.
    void addEdge(
        const VectorMax3d& edge_vertex0,
        const VectorMax3d& edge_vertex1,
        const long index,
        const double inflation_radius = 0.0)
    {
        addEdge(
            edge_vertex0, edge_vertex1, edge_vertex0, edge_vertex1, index,
            inflation_radius);
    }

    /// @brief Add all edges as AABBs containing the time swept quad.
    void addEdges(
        const Eigen::MatrixXd& vertices_t0,
        const Eigen::MatrixXd& vertices_t1,
        const Eigen::MatrixXi& edges,
        const double inflation_radius = 0.0);

    /// @brief Add all edges as AABBs.
    void addEdges(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const double inflation_radius = 0.0)
    {
        addEdges(vertices, vertices, edges, inflation_radius);
    }

    /// @brief Add an edge as a AABB containing the time swept prism.
    void addFace(
        const VectorMax3d& face_vertex0_t0,
        const VectorMax3d& face_vertex1_t0,
        const VectorMax3d& face_vertex2_t0,
        const VectorMax3d& face_vertex0_t1,
        const VectorMax3d& face_vertex1_t1,
        const VectorMax3d& face_vertex2_t1,
        const long index,
        const double inflation_radius = 0.0);

    /// @brief Add an edge as a AABB.
    void addFace(
        const VectorMax3d& face_vertex0,
        const VectorMax3d& face_vertex1,
        const VectorMax3d& face_vertex2,
        const long index,
        const double inflation_radius = 0.0)
    {
        addFace(
            face_vertex0, face_vertex1, face_vertex2, face_vertex0,
            face_vertex1, face_vertex2, index, inflation_radius);
    }

    /// @brief Add all edges as AABBs containing the time swept prism.
    void addFaces(
        const Eigen::MatrixXd& vertices_t0,
        const Eigen::MatrixXd& vertices_t1,
        const Eigen::MatrixXi& faces,
        const double inflation_radius = 0.0);

    /// @brief Add all edges as AABBs.
    void addFaces(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& faces,
        const double inflation_radius = 0.0)
    {
        addFaces(vertices, vertices, faces, inflation_radius);
    }

    /// @brief Compute the candidate edge-vertex candidate collisisons.
    void getVertexEdgePairs(
        const Eigen::MatrixXi& edges,
        std::vector<EdgeVertexCandidate>& ev_candidates,
        const std::function<bool(size_t, size_t)>& can_collide =
            [](size_t, size_t) { return true; });

    /// @brief Compute the candidate edge-edge candidate collisions.
    void getEdgeEdgePairs(
        const Eigen::MatrixXi& edges,
        std::vector<EdgeEdgeCandidate>& ee_candidates,
        const std::function<bool(size_t, size_t)>& can_collide =
            [](size_t, size_t) { return true; });

    /// @brief Compute the candidate edge-face candidate intersections.
    void getEdgeFacePairs(
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        std::vector<EdgeFaceCandidate>& ef_candidates,
        const std::function<bool(size_t, size_t)>& can_collide =
            [](size_t, size_t) { return true; });

    /// @brief Compute the candidate edge-edge candidate collisions.
    void getFaceVertexPairs(
        const Eigen::MatrixXi& faces,
        std::vector<FaceVertexCandidate>& fv_candidates,
        const std::function<bool(size_t, size_t)>& can_collide =
            [](size_t, size_t) { return true; });

    /// @brief Clear the hash grid.
    inline void clear()
    {
        m_vertexItems.clear();
        m_edgeItems.clear();
        m_faceItems.clear();
    }

protected:
    /// @brief Add an AABB of the extents to the hash grid.
    void addElement(
        const AABB& aabb, const int id, std::vector<HashItem>& items) const;

    /// @brief Add a vertex as a AABB containing the time swept edge.
    void addVertex(
        const VectorMax3d& vertex_t0,
        const VectorMax3d& vertex_t1,
        const long index,
        std::vector<HashItem>& vertex_items,
        const double inflation_radius = 0.0) const;

    /// @brief Add an edge as a AABB containing the time swept quad.
    void addEdge(
        const VectorMax3d& edge_vertex0_t0,
        const VectorMax3d& edge_vertex1_t0,
        const VectorMax3d& edge_vertex0_t1,
        const VectorMax3d& edge_vertex1_t1,
        const long index,
        std::vector<HashItem>& edge_items,
        const double inflation_radius = 0.0) const;

    /// @brief Add an edge as a AABB containing the time swept quad.
    void addFace(
        const VectorMax3d& face_vertex0_t0,
        const VectorMax3d& face_vertex1_t0,
        const VectorMax3d& face_vertex2_t0,
        const VectorMax3d& face_vertex0_t1,
        const VectorMax3d& face_vertex1_t1,
        const VectorMax3d& face_vertex2_t1,
        const long index,
        std::vector<HashItem>& face_items,
        const double inflation_radius = 0.0) const;

    /// @brief Create the hash of a cell location.
    inline long hash(int x, int y, int z) const
    {
        assert(x >= 0 && y >= 0 && z >= 0);
        assert(
            x < m_gridSize[0] && y < m_gridSize[1]
            && (m_gridSize.size() == 2 || z < m_gridSize[2]));
        return (z * m_gridSize[1] + y) * m_gridSize[0] + x;
    }

protected:
    double m_cellSize;
    ArrayMax3i m_gridSize;
    ArrayMax3d m_domainMin;
    ArrayMax3d m_domainMax;

    std::vector<HashItem> m_vertexItems;
    std::vector<HashItem> m_edgeItems;
    std::vector<HashItem> m_faceItems;
};

} // namespace ipc
