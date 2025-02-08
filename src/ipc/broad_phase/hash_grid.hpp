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
    HashGrid() = default;

    /// @brief Get the name of the broad phase method.
    /// @return The name of the broad phase method.
    std::string name() const override { return "HashGrid"; }

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

    /// @brief Find the candidate face-face collisions.
    /// @param[out] candidates The candidate face-face collisions.
    void detect_face_face_candidates(
        std::vector<FaceFaceCandidate>& candidates) const override;

    double cell_size() const { return m_cell_size; }
    const ArrayMax3i& grid_size() const { return m_grid_size; }
    const ArrayMax3d& domain_min() const { return m_domain_min; }
    const ArrayMax3d& domain_max() const { return m_domain_max; }

protected:
    void resize(
        const ArrayMax3d& domain_min,
        const ArrayMax3d& domain_max,
        double cell_size);

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
            x < grid_size()[0] && y < grid_size()[1]
            && (grid_size().size() == 2 || z < grid_size()[2]));
        return (z * grid_size()[1] + y) * grid_size()[0] + x;
    }

private:
    /// @brief Find the candidate collisions between two sets of items.
    /// @tparam Candidate The type of collision candidate.
    /// @param[in] items0 First set of items.
    /// @param[in] items1 Second set of items.
    /// @param[in] boxes0 First set's boxes.
    /// @param[in] boxes1 Second set's boxes.
    /// @param[in] can_collide Function to determine if two items can collide.
    /// @param[out] candidates The candidate collisions.
    template <typename Candidate>
    void detect_candidates(
        const std::vector<HashItem>& items0,
        const std::vector<HashItem>& items1,
        const std::vector<AABB>& boxes0,
        const std::vector<AABB>& boxes1,
        const std::function<bool(size_t, size_t)>& can_collide,
        std::vector<Candidate>& candidates) const;

    /// @brief Find the candidate collisions among a set of items.
    /// @tparam Candidate The type of collision candidate.
    /// @param[in] items The set of items.
    /// @param[in] boxes The items' boxes.
    /// @param[in] can_collide Function to determine if two items can collide.
    /// @param[out] candidates The candidate collisions.
    template <typename Candidate>
    void detect_candidates(
        const std::vector<HashItem>& items,
        const std::vector<AABB>& boxes,
        const std::function<bool(size_t, size_t)>& can_collide,
        std::vector<Candidate>& candidates) const;

protected:
    double m_cell_size;
    ArrayMax3i m_grid_size;
    ArrayMax3d m_domain_min;
    ArrayMax3d m_domain_max;

    std::vector<HashItem> vertex_items;
    std::vector<HashItem> edge_items;
    std::vector<HashItem> face_items;
};

} // namespace ipc
