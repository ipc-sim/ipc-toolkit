#include <ipc/spatial_hash/hash_grid.hpp>

#ifdef IPC_TOOLKIT_SPATIAL_HASH_USE_TBB
#include <tbb/parallel_for.h>
#endif
#include <tbb/parallel_sort.h> // Still use this even if TBB is disabled

#ifdef IPC_TOOLKIT_WITH_LOGGER
#include <ipc/utils/logger.hpp>
#endif
namespace ipc {

bool AABB::are_overlapping(const AABB& a, const AABB& b)
{
    // https://bit.ly/2ZP3tW4
    assert(a.dim == b.dim);
    return (abs(a.center.x() - b.center.x())
            <= (a.half_extent.x() + b.half_extent.x()))
        && (abs(a.center.y() - b.center.y())
            <= (a.half_extent.y() + b.half_extent.y()))
        && (a.dim == 2
            || abs(a.center.z() - b.center.z())
                <= (a.half_extent.z() + b.half_extent.z()));
};

void HashGrid::resize(
    Eigen::VectorX3d min, Eigen::VectorX3d max, double cellSize)
{
    clear();
    assert(cellSize != 0.0);
    m_cellSize = cellSize;
    m_domainMin = min;
    m_domainMax = max;
    m_gridSize = ((max - min) / m_cellSize).array().ceil().cast<int>().max(1);
#ifdef IPC_TOOLKIT_WITH_LOGGER
    logger().debug(
        "hash-grid resized with a size of {:d}x{:d}x{:d}", m_gridSize[0],
        m_gridSize[1], m_gridSize.size() == 3 ? m_gridSize[2] : 1);
#endif
}

/// @brief Compute an AABB around a given 2D mesh.
void calculate_mesh_extents(
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    Eigen::VectorX3d& lower_bound,
    Eigen::VectorX3d& upper_bound)
{
    int dim = vertices_t0.cols();
    Eigen::MatrixXd points(vertices_t0.rows() + vertices_t1.rows(), dim);
    points.topRows(vertices_t0.rows()) = vertices_t0;
    points.bottomRows(vertices_t1.rows()) = vertices_t1;

    lower_bound = points.colwise().minCoeff();
    upper_bound = points.colwise().maxCoeff();
}

/// @brief Compute the average edge length of a mesh.
double average_edge_length(
    const Eigen::MatrixXd& V_t0,
    const Eigen::MatrixXd& V_t1,
    const Eigen::MatrixXi& E)
{
    double avg = 0;
    for (unsigned i = 0; i < E.rows(); ++i) {
        avg += (V_t0.row(E(i, 0)) - V_t0.row(E(i, 1))).norm();
        avg += (V_t1.row(E(i, 0)) - V_t1.row(E(i, 1))).norm();
    }
    return avg / (2 * E.rows());
}

/// @brief Compute the average displacement length.
double average_displacement_length(const Eigen::MatrixXd& displacements)
{
    return displacements.rowwise().norm().sum() / displacements.rows();
}

void HashGrid::resize(
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const Eigen::MatrixXi& edges,
    const double inflation_radius)
{
    Eigen::VectorX3d mesh_min, mesh_max;
    calculate_mesh_extents(vertices_t0, vertices_t1, mesh_min, mesh_max);
    double edge_len = average_edge_length(vertices_t0, vertices_t1, edges);
    double disp_len = average_displacement_length(vertices_t1 - vertices_t0);
    double cell_size = 2 * std::max(edge_len, disp_len) + inflation_radius;
    this->resize(
        mesh_min.array() - inflation_radius,
        mesh_max.array() + inflation_radius, cell_size);
}

/// @brief Compute a AABB for a vertex moving through time (i.e. temporal edge).
void calculate_vertex_extents(
    const Eigen::VectorX3d& vertex_t0,
    const Eigen::VectorX3d& vertex_t1,
    Eigen::VectorX3d& lower_bound,
    Eigen::VectorX3d& upper_bound)
{
    Eigen::MatrixXd points(2, vertex_t0.size());
    points.row(0) = vertex_t0;
    points.row(1) = vertex_t1;

    lower_bound = points.colwise().minCoeff();
    upper_bound = points.colwise().maxCoeff();
}

void HashGrid::addVertex(
    const Eigen::VectorX3d& vertex_t0,
    const Eigen::VectorX3d& vertex_t1,
    const long index,
    const double inflation_radius)
{
    Eigen::VectorX3d lower_bound, upper_bound;
    calculate_vertex_extents(vertex_t0, vertex_t1, lower_bound, upper_bound);
    this->addElement(
        AABB(
            lower_bound.array() - inflation_radius,
            upper_bound.array() + inflation_radius),
        index, m_vertexItems); // Vertices have a negative id
}

void HashGrid::addVertices(
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const double inflation_radius)
{
    assert(vertices_t0.rows() == vertices_t1.rows());
#ifdef IPC_TOOLKIT_SPATIAL_HASH_USE_TBB
    tbb::parallel_for(0l, (long)(vertices_t0.rows()), [&](long i) {
#else
    for (long i = 0; i < vertices_t0.rows(); i++) {
#endif
        addVertex(vertices_t0.row(i), vertices_t1.row(i), i, inflation_radius);
#ifdef IPC_TOOLKIT_SPATIAL_HASH_USE_TBB
    });
#else
    }
#endif
}

/// @brief Compute a AABB for an edge moving through time (i.e. temporal quad).
void calculate_edge_extents(
    const Eigen::VectorX3d& edge_vertex0_t0,
    const Eigen::VectorX3d& edge_vertex1_t0,
    const Eigen::VectorX3d& edge_vertex0_t1,
    const Eigen::VectorX3d& edge_vertex1_t1,
    Eigen::VectorX3d& lower_bound,
    Eigen::VectorX3d& upper_bound)
{
    Eigen::MatrixXd points(4, edge_vertex0_t0.size());
    points.row(0) = edge_vertex0_t0;
    points.row(1) = edge_vertex1_t0;
    points.row(2) = edge_vertex0_t1;
    points.row(3) = edge_vertex1_t1;

    lower_bound = points.colwise().minCoeff();
    upper_bound = points.colwise().maxCoeff();
}

void HashGrid::addEdge(
    const Eigen::VectorX3d& edge_vertex0_t0,
    const Eigen::VectorX3d& edge_vertex1_t0,
    const Eigen::VectorX3d& edge_vertex0_t1,
    const Eigen::VectorX3d& edge_vertex1_t1,
    const long index,
    const double inflation_radius)
{
    Eigen::VectorX3d lower_bound, upper_bound;
    calculate_edge_extents(
        edge_vertex0_t0, edge_vertex1_t0, edge_vertex0_t1, edge_vertex1_t1,
        lower_bound, upper_bound);
    this->addElement(
        AABB(
            lower_bound.array() - inflation_radius,
            upper_bound.array() + inflation_radius),
        index, m_edgeItems); // Edges have a positive id
}

void HashGrid::addEdges(
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const Eigen::MatrixXi& edges,
    const double inflation_radius)
{
    assert(vertices_t0.rows() == vertices_t1.rows());
#ifdef IPC_TOOLKIT_SPATIAL_HASH_USE_TBB
    tbb::parallel_for(0l, (long)(edges.rows()), [&](long i) {
#else
    for (long i = 0; i < edges.rows(); i++) {
#endif
        addEdge(
            vertices_t0.row(edges(i, 0)), vertices_t0.row(edges(i, 1)),
            vertices_t1.row(edges(i, 0)), vertices_t1.row(edges(i, 1)), i,
            inflation_radius);
#ifdef IPC_TOOLKIT_SPATIAL_HASH_USE_TBB
    });
#else
    }
#endif
}

/// @brief Compute a AABB for an edge moving through time (i.e. temporal quad).
void calculate_face_extents(
    const Eigen::VectorX3d& face_vertex0_t0,
    const Eigen::VectorX3d& face_vertex1_t0,
    const Eigen::VectorX3d& face_vertex2_t0,
    const Eigen::VectorX3d& face_vertex0_t1,
    const Eigen::VectorX3d& face_vertex1_t1,
    const Eigen::VectorX3d& face_vertex2_t1,
    Eigen::VectorX3d& lower_bound,
    Eigen::VectorX3d& upper_bound)
{
    Eigen::MatrixXd points(6, face_vertex0_t0.size());
    points.row(0) = face_vertex0_t0;
    points.row(1) = face_vertex1_t0;
    points.row(2) = face_vertex2_t0;
    points.row(3) = face_vertex0_t1;
    points.row(4) = face_vertex1_t1;
    points.row(5) = face_vertex2_t1;

    lower_bound = points.colwise().minCoeff();
    upper_bound = points.colwise().maxCoeff();
}

void HashGrid::addFace(
    const Eigen::VectorX3d& face_vertex0_t0,
    const Eigen::VectorX3d& face_vertex1_t0,
    const Eigen::VectorX3d& face_vertex2_t0,
    const Eigen::VectorX3d& face_vertex0_t1,
    const Eigen::VectorX3d& face_vertex1_t1,
    const Eigen::VectorX3d& face_vertex2_t1,
    const long index,
    const double inflation_radius)
{
    Eigen::VectorX3d lower_bound, upper_bound;
    calculate_face_extents(
        face_vertex0_t0, face_vertex1_t0, face_vertex2_t0, //
        face_vertex0_t1, face_vertex1_t1, face_vertex2_t1, //
        lower_bound, upper_bound);
    this->addElement(
        AABB(
            lower_bound.array() - inflation_radius,
            upper_bound.array() + inflation_radius),
        index, m_faceItems); // Faces have a positive id
}

void HashGrid::addFaces(
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const Eigen::MatrixXi& faces,
    const double inflation_radius)
{
    assert(vertices_t0.rows() == vertices_t1.rows());
#ifdef IPC_TOOLKIT_SPATIAL_HASH_USE_TBB
    tbb::parallel_for(0l, (long)(faces.rows()), [&](long i) {
#else
    for (long i = 0; i < faces.rows(); i++) {
#endif
        addFace(
            vertices_t0.row(faces(i, 0)), vertices_t0.row(faces(i, 1)),
            vertices_t0.row(faces(i, 2)), vertices_t1.row(faces(i, 0)),
            vertices_t1.row(faces(i, 1)), vertices_t1.row(faces(i, 2)), i,
            inflation_radius);

#ifdef IPC_TOOLKIT_SPATIAL_HASH_USE_TBB
    });
#else
    }
#endif
}

void HashGrid::addElement(const AABB& aabb, const int id, HashItems& items)
{
    Eigen::VectorX3<int> int_min =
        ((aabb.getMin() - m_domainMin) / m_cellSize).cast<int>();
    // We can round down to -1, but not less
    assert((int_min.array() >= -1).all());
    assert((int_min.array() <= m_gridSize.array()).all());
    int_min = int_min.array().max(0).min(m_gridSize.array() - 1);

    Eigen::VectorX3<int> int_max =
        ((aabb.getMax() - m_domainMin) / m_cellSize).cast<int>();
    assert((int_max.array() >= -1).all());
    assert((int_max.array() <= m_gridSize.array()).all());
    int_max = int_max.array().max(0).min(m_gridSize.array() - 1);
    assert((int_min.array() <= int_max.array()).all());

    int min_z = int_min.size() == 3 ? int_min.z() : 0;
    int max_z = int_max.size() == 3 ? int_max.z() : 0;
    for (int x = int_min.x(); x <= int_max.x(); ++x) {
        for (int y = int_min.y(); y <= int_max.y(); ++y) {
            for (int z = min_z; z <= max_z; ++z) {
                items.emplace_back(hash(x, y, z), id, aabb);
            }
        }
    }
}

template <typename T>
void getPairs(
    const std::function<bool(int, int)>& is_endpoint,
    const std::function<bool(int, int)>& is_same_group,
    HashItems& items0,
    HashItems& items1,
    T& candidates)
{
    // Sorted all they (key,value) pairs, where key is the hash key, and value
    // is the element index
    tbb::parallel_sort(items0.begin(), items0.end());
    tbb::parallel_sort(items1.begin(), items1.end());

    // Entries with the same key means they share a cell (that cell index
    // hashes to the same key) and should be flagged for low-level intersection
    // testing. So we loop over the entire sorted set of (key,value) pairs
    // creating Candidate entries for vertex-edge pairs with the same key
    int i = 0, j_start = 0;
    while (i < items0.size() && j_start < items1.size()) {
        const HashItem& item0 = items0[i];

        int j = j_start;
        while (j < items1.size()) {
            const HashItem& item1 = items1[j];

            if (item0.key == item1.key) {
                if (!is_endpoint(item0.id, item1.id)
                    && !is_same_group(item0.id, item1.id)
                    && AABB::are_overlapping(item0.aabb, item1.aabb)) {
                    candidates.emplace_back(item0.id, item1.id);
                }
            } else {
                break;
            }
            j++;
        }

        if (i == items0.size() - 1 || item0.key != items0[i + 1].key) {
            j_start = j;
        }
        i++;
    }

    // Remove the duplicate candidates
    tbb::parallel_sort(candidates.begin(), candidates.end());
    auto new_end = std::unique(candidates.begin(), candidates.end());
    candidates.erase(new_end, candidates.end());
}

template <typename T>
void getPairs(
    const std::function<bool(int, int)>& is_endpoint,
    const std::function<bool(int, int)>& is_same_group,
    HashItems& items,
    T& candidates)
{
    // Sorted all they (key,value) pairs, where key is the hash key, and value
    // is the element index
    tbb::parallel_sort(items.begin(), items.end());

    // Entries with the same key means they share a cell (that cell index
    // hashes to the same key) and should be flagged for low-level intersection
    // testing. So we loop over the entire sorted set of (key,value) pairs
    // creating Candidate entries for vertex-edge pairs with the same key
    for (int i = 0; i < items.size(); i++) {
        const HashItem& item0 = items[i];
        for (int j = i + 1; j < items.size(); j++) {
            const HashItem& item1 = items[j];
            if (item0.key == item1.key) {
                if (!is_endpoint(item0.id, item1.id)
                    && !is_same_group(item0.id, item1.id)
                    && AABB::are_overlapping(item0.aabb, item1.aabb)) {
                    candidates.emplace_back(item0.id, item1.id);
                }
            } else {
                break; // This avoids a brute force comparison
            }
        }
    }

    // Remove the duplicate candidates
    tbb::parallel_sort(candidates.begin(), candidates.end());
    auto new_end = std::unique(candidates.begin(), candidates.end());
    candidates.erase(new_end, candidates.end());
}

void HashGrid::getVertexEdgePairs(
    const Eigen::MatrixXi& edges,
    const Eigen::VectorXi& group_ids,
    std::vector<EdgeVertexCandidate>& ev_candidates)
{
    auto is_endpoint = [&](int ei, int vi) {
        return edges(ei, 0) == vi || edges(ei, 1) == vi;
    };

    bool check_groups = group_ids.size() > 0;
    auto is_same_group = [&](int ei, int vi) {
        return check_groups
            && (group_ids(vi) == group_ids(edges(ei, 0))
                || group_ids(vi) == group_ids(edges(ei, 1)));
    };

    getPairs(
        is_endpoint, is_same_group, m_edgeItems, m_vertexItems, ev_candidates);
}

void HashGrid::getEdgeEdgePairs(
    const Eigen::MatrixXi& edges,
    const Eigen::VectorXi& group_ids,
    std::vector<EdgeEdgeCandidate>& ee_candidates)
{
    auto is_endpoint = [&](int ei, int ej) {
        return edges(ei, 0) == edges(ej, 0) || edges(ei, 0) == edges(ej, 1)
            || edges(ei, 1) == edges(ej, 0) || edges(ei, 1) == edges(ej, 1);
    };

    bool check_groups = group_ids.size() > 0;
    auto is_same_group = [&](int ei, int ej) {
        return check_groups
            && (group_ids(edges(ei, 0)) == group_ids(edges(ej, 0))
                || group_ids(edges(ei, 0)) == group_ids(edges(ej, 1))
                || group_ids(edges(ei, 1)) == group_ids(edges(ej, 0))
                || group_ids(edges(ei, 1)) == group_ids(edges(ej, 1)));
    };

    getPairs(is_endpoint, is_same_group, m_edgeItems, ee_candidates);
}

void HashGrid::getEdgeFacePairs(
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const Eigen::VectorXi& group_ids,
    std::vector<EdgeFaceCandidate>& ef_candidates)
{
    auto is_endpoint = [&](int ei, int fi) {
        // Check if the edge and face have a common end-point
        return edges(ei, 0) == faces(fi, 0) || edges(ei, 0) == faces(fi, 1)
            || edges(ei, 0) == faces(fi, 2) || edges(ei, 1) == faces(fi, 0)
            || edges(ei, 1) == faces(fi, 1) || edges(ei, 1) == faces(fi, 2);
    };

    bool check_groups = group_ids.size() > 0;
    auto is_same_group = [&](int ei, int fi) {
        return check_groups
            && (group_ids(edges(ei, 0)) == group_ids(faces(fi, 0))
                || group_ids(edges(ei, 0)) == group_ids(faces(fi, 1))
                || group_ids(edges(ei, 0)) == group_ids(faces(fi, 2))
                || group_ids(edges(ei, 1)) == group_ids(faces(fi, 0))
                || group_ids(edges(ei, 1)) == group_ids(faces(fi, 1))
                || group_ids(edges(ei, 1)) == group_ids(faces(fi, 2)));
    };

    getPairs(
        is_endpoint, is_same_group, m_edgeItems, m_faceItems, ef_candidates);
}

void HashGrid::getFaceVertexPairs(
    const Eigen::MatrixXi& faces,
    const Eigen::VectorXi& group_ids,
    std::vector<FaceVertexCandidate>& fv_candidates)
{
    auto is_endpoint = [&](int fi, int vi) {
        return vi == faces(fi, 0) || vi == faces(fi, 1) || vi == faces(fi, 2);
    };

    bool check_groups = group_ids.size() > 0;
    auto is_same_group = [&](int fi, int vi) {
        return check_groups
            && (group_ids(vi) == group_ids(faces(fi, 0))
                || group_ids(vi) == group_ids(faces(fi, 1))
                || group_ids(vi) == group_ids(faces(fi, 2)));
    };

    getPairs(
        is_endpoint, is_same_group, m_faceItems, m_vertexItems, fv_candidates);
}

} // namespace ipc
