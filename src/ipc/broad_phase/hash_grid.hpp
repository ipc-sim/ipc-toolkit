#pragma once

#include <ipc/broad_phase/broad_phase.hpp>

namespace ipc {

/// @brief An entry into the hash grid as a (key, value) pair.
struct HashItem {
    /// @brief The key of the item.
    long key;
    /// @brief The value of the item.
    long id;

    /// @brief Construct a hash item as a (key, value) pair.
    HashItem(int _key, int _id) : key(_key), id(_id) { }

    /// @brief Compare HashItems by their keys for sorting.
    bool operator<(const HashItem& other) const
    {
        if (key == other.key) {
            return id < other.id;
        }
        return key < other.key;
    }
};

class HashGrid : public BroadPhase {
public:
    /// @brief Build the broad phase for static collision detection.
    /// @param vertices Vertex positions
    /// @param edges Collision mesh edges
    /// @param faces Collision mesh faces
    /// @param inflation_radius Radius of inflation around all elements.
    void build(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        double inflation_radius = 0) override;

    /// @brief Build the broad phase for continuous collision detection.
    /// @param vertices_t0 Starting vertices of the vertices.
    /// @param vertices_t1 Ending vertices of the vertices.
    /// @param edges Collision mesh edges
    /// @param faces Collision mesh faces
    /// @param inflation_radius Radius of inflation around all elements.
    void build(
        const Eigen::MatrixXd& vertices_t0,
        const Eigen::MatrixXd& vertices_t1,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        double inflation_radius = 0) override;

    /// @brief Clear the hash grid.
    void clear() override
    {
        BroadPhase::clear();
        vertex_items.clear();
        edge_items.clear();
        face_items.clear();
    }

    /// @brief Find the candidate vertex-vertex collisions.
    void detect_vertex_vertex_candidates(
        std::vector<VertexVertexCandidate>& candidates) const override;

    /// @brief Find the candidate edge-vertex collisions.
    /// @param[out] candidates The candidate edge-vertex collisions.
    void detect_edge_vertex_candidates(
        std::vector<EdgeVertexCandidate>& candidates) const override;

    /// @brief Find the candidate edge-edge collisions.
    /// @param[out] candidates The candidate edge-edge collisions.
    void detect_edge_edge_candidates(
        std::vector<EdgeEdgeCandidate>& candidates) const override;

    /// @brief Find the candidate face-vertex collisions.
    /// @param[out] candidates The candidate face-vertex collisions.
    void detect_face_vertex_candidates(
        std::vector<FaceVertexCandidate>& candidates) const override;

    /// @brief Find the candidate edge-face intersections.
    /// @param[out] candidates The candidate edge-face intersections.
    void detect_edge_face_candidates(
        std::vector<EdgeFaceCandidate>& candidates) const override;

    double cellSize() const { return m_cellSize; }
    const ArrayMax3i& gridSize() const { return m_gridSize; }
    const ArrayMax3d& domainMin() const { return m_domainMin; }
    const ArrayMax3d& domainMax() const { return m_domainMax; }

protected:
    void resize(const ArrayMax3d& min, const ArrayMax3d& max, double cellSize);

    void insert_boxes();

    void insert_boxes(
        const std::vector<AABB>& boxes, std::vector<HashItem>& items) const;

    /// @brief Add an AABB of the extents to the hash grid.
    void insert_box(
        const AABB& aabb, const long id, std::vector<HashItem>& items) const;

    /// @brief Create the hash of a cell location.
    inline long hash(int x, int y, int z) const
    {
        assert(x >= 0 && y >= 0 && z >= 0);
        assert(
            x < m_gridSize[0] && y < m_gridSize[1]
            && (m_gridSize.size() == 2 || z < m_gridSize[2]));
        return (z * m_gridSize[1] + y) * m_gridSize[0] + x;
    }

private:
    template <typename Candidate>
    void detect_candidates(
        const std::vector<HashItem>& items0,
        const std::vector<HashItem>& items1,
        const std::vector<AABB>& boxes0,
        const std::vector<AABB>& boxes1,
        const std::function<bool(size_t, size_t)>& can_collide,
        std::vector<Candidate>& candidates) const;

    template <typename Candidate>
    void detect_candidates(
        const std::vector<HashItem>& items,
        const std::vector<AABB>& boxes,
        const std::function<bool(size_t, size_t)>& can_collide,
        std::vector<Candidate>& candidates) const;

protected:
    double m_cellSize;
    ArrayMax3i m_gridSize;
    ArrayMax3d m_domainMin;
    ArrayMax3d m_domainMax;

    std::vector<HashItem> vertex_items;
    std::vector<HashItem> edge_items;
    std::vector<HashItem> face_items;
};

} // namespace ipc
