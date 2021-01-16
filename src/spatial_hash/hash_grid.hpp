#pragma once

#include <vector>

#include <Eigen/Core>

#include <ipc/spatial_hash/collision_candidate.hpp>
#include <ipc/utils/eigen_ext.hpp>

// WARNING: In limited testing the way TBB is used in the hash grid was shown to
// be detremental to performance, so we currently recommend leaving this line
// commented out.
// TODO: Implementing a more parallel or integrating the original spatial hash
// is an important goal in order to increase performance.
// #define IPC_TOOLKIT_SPATIAL_HASH_USE_TBB
#ifdef IPC_TOOLKIT_SPATIAL_HASH_USE_TBB
#include <tbb/concurrent_vector.h>
#endif

namespace ipc {

/// @brief Axis aligned bounding-box of some type
class AABB {
public:
    AABB() {}

    AABB(const Eigen::VectorX3d& min, const Eigen::VectorX3d& max)
        : min(min)
        , max(max)
    {
        assert((min.array() <= max.array()).all());
        half_extent = (max - min) / 2;
        center = min + half_extent;
        dim = min.size();
        assert(max.size() == dim);
    }

    AABB(const AABB& aabb1, const AABB& aabb2)
        : AABB(
              aabb1.min.array().min(aabb2.min.array()),
              aabb1.max.array().max(aabb2.max.array()))
    {
    }

    AABB(const AABB& aabb1, const AABB& aabb2, const AABB& aabb3)
        : AABB(AABB(aabb1, aabb2), aabb3)
    {
    }

    static bool are_overlapping(const AABB& a, const AABB& b);

    inline const Eigen::VectorX3d& getMin() const { return min; }
    inline const Eigen::VectorX3d& getMax() const { return max; }
    inline const Eigen::VectorX3d& getHalfExtent() const { return half_extent; }
    inline const Eigen::VectorX3d& getCenter() const { return center; }

protected:
    Eigen::VectorX3d min;
    Eigen::VectorX3d max;
    Eigen::VectorX3d half_extent;
    Eigen::VectorX3d center;
    int dim;
};

/// @brief An entry into the hash grid as a (key, value) pair.
class HashItem {
public:
    int key;   /// @brief The key of the item.
    int id;    /// @brief The value of the item.
    AABB aabb; /// @brief The axis-aligned bounding box of the element

    /// @brief Construct a hash item as a (key, value) pair.
    HashItem(int key, int id, const AABB aabb)
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

#ifdef IPC_TOOLKIT_SPATIAL_HASH_USE_TBB
// TODO: This may be less efficient than a std::vector
typedef tbb::concurrent_vector<HashItem> HashItems;
#else
typedef std::vector<HashItem> HashItems;
#endif

class HashGrid {
public:
    void resize(Eigen::VectorX3d min, Eigen::VectorX3d max, double cellSize);

    void resize(
        const Eigen::MatrixXd& vertices_t0,
        const Eigen::MatrixXd& vertices_t1,
        const Eigen::MatrixXi& edges,
        const double inflation_radius = 0.0);

    /// @brief Add a vertex as a AABB containing the time swept edge.
    void addVertex(
        const Eigen::VectorX3d& vertex_t0,
        const Eigen::VectorX3d& vertex_t1,
        const long index,
        const double inflation_radius = 0.0);

    /// @brief Add all vertices as AABBs containing the time swept edge.
    void addVertices(
        const Eigen::MatrixXd& vertices_t0,
        const Eigen::MatrixXd& vertices_t1,
        const double inflation_radius = 0.0);

    /// @brief Add an edge as a AABB containing the time swept quad.
    void addEdge(
        const Eigen::VectorX3d& edge_vertex0_t0,
        const Eigen::VectorX3d& edge_vertex1_t0,
        const Eigen::VectorX3d& edge_vertex0_t1,
        const Eigen::VectorX3d& edge_vertex1_t1,
        const long index,
        const double inflation_radius = 0.0);

    /// @brief Add all edges as AABBs containing the time swept quad.
    void addEdges(
        const Eigen::MatrixXd& vertices_t0,
        const Eigen::MatrixXd& vertices_t1,
        const Eigen::MatrixXi& edges,
        const double inflation_radius = 0.0);

    /// @brief Add an edge as a AABB containing the time swept quad.
    void addFace(
        const Eigen::VectorX3d& face_vertex0_t0,
        const Eigen::VectorX3d& face_vertex1_t0,
        const Eigen::VectorX3d& face_vertex2_t0,
        const Eigen::VectorX3d& face_vertex0_t1,
        const Eigen::VectorX3d& face_vertex1_t1,
        const Eigen::VectorX3d& face_vertex2_t1,
        const long index,
        const double inflation_radius = 0.0);

    /// @brief Add all edges as AABBs containing the time swept quad.
    void addFaces(
        const Eigen::MatrixXd& vertices_t0,
        const Eigen::MatrixXd& vertices_t1,
        const Eigen::MatrixXi& faces,
        const double inflation_radius = 0.0);

    /// @brief Compute the candidate edge-vertex candidate collisisons.
    void getVertexEdgePairs(
        const Eigen::MatrixXi& edges,
        const Eigen::VectorXi& group_ids,
        std::vector<EdgeVertexCandidate>& ev_candidates);

    /// @brief Compute the candidate edge-edge candidate collisions.
    void getEdgeEdgePairs(
        const Eigen::MatrixXi& edges,
        const Eigen::VectorXi& group_ids,
        std::vector<EdgeEdgeCandidate>& ee_candidates);

    /// @brief Compute the candidate edge-face candidate intersections.
    void getEdgeFacePairs(
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        const Eigen::VectorXi& group_ids,
        std::vector<EdgeFaceCandidate>& ef_candidates);

    /// @brief Compute the candidate edge-edge candidate collisions.
    void getFaceVertexPairs(
        const Eigen::MatrixXi& faces,
        const Eigen::VectorXi& group_ids,
        std::vector<FaceVertexCandidate>& fv_candidates);

    /// @brief Clear the hash grid.
    inline void clear()
    {
        m_vertexItems.clear();
        m_edgeItems.clear();
        m_faceItems.clear();
    }

protected:
    /// @brief Add an AABB of the extents to the hash grid.
    void addElement(const AABB& aabb, const int id, HashItems& items);

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
    Eigen::VectorX3i m_gridSize;
    Eigen::VectorX3d m_domainMin;
    Eigen::VectorX3d m_domainMax;

    HashItems m_vertexItems;
    HashItems m_edgeItems;
    HashItems m_faceItems;
};

} // namespace ipc
