#include "collision_mesh.hpp"

#include <ipc/distance/point_edge.hpp>

#include <ipc/utils/unordered_map_and_set.hpp>
#include <ipc/utils/logger.hpp>
#include <ipc/utils/eigen_ext.hpp>
#include <ipc/utils/math.hpp>
#include <ipc/utils/local_to_global.hpp>
#include <ipc/utils/area_gradient.hpp>

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

namespace ipc {

CollisionMesh::CollisionMesh(
    const Eigen::MatrixXd& rest_positions,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const Eigen::SparseMatrix<double>& displacement_map)
    : CollisionMesh(
        std::vector<bool>(rest_positions.rows(), true),
        rest_positions,
        edges,
        faces,
        displacement_map)
{
}

CollisionMesh::CollisionMesh(
    const std::vector<bool>& include_vertex,
    const Eigen::MatrixXd& full_rest_positions,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const Eigen::SparseMatrix<double>& displacement_map)
    : m_full_rest_positions(full_rest_positions)
    , m_edges(edges)
    , m_faces(faces)
{
    assert(include_vertex.size() == full_rest_positions.rows());
    const bool include_all_vertices = std::all_of(
        include_vertex.begin(), include_vertex.end(), [](bool b) { return b; });

    if (include_all_vertices) {
        // set full ↔ reduced ≡ identity
        m_full_vertex_to_vertex.setLinSpaced(
            full_rest_positions.rows(), 0, full_rest_positions.rows() - 1);
        m_vertex_to_full_vertex = m_full_vertex_to_vertex;
    } else {
        m_full_vertex_to_vertex.setConstant(full_rest_positions.rows(), -1);
        std::vector<int> dynamic_vertex_to_full_vertex;
        for (size_t i = 0; i < full_rest_positions.rows(); i++) {
            if (include_vertex[i]) {
                m_full_vertex_to_vertex[i] =
                    dynamic_vertex_to_full_vertex.size();
                dynamic_vertex_to_full_vertex.push_back(i);
            }
        }
        m_vertex_to_full_vertex = Eigen::Map<Eigen::VectorXi>(
            dynamic_vertex_to_full_vertex.data(),
            dynamic_vertex_to_full_vertex.size());
    }

    // ========================================================================

    const int dim = full_rest_positions.cols();

    // Selection matrix S ∈ ℝ^{collision×full}
    init_selection_matrices(dim);

    if (displacement_map.size() == 0) {
        m_displacement_map = m_select_vertices;
        m_displacement_dof_map = m_select_dof;
    } else {
        assert(displacement_map.rows() == full_num_vertices());
        // assert(displacement_map.cols() == full_num_vertices());

        m_displacement_map = m_select_vertices * displacement_map;
        m_displacement_map.makeCompressed();

        m_displacement_dof_map =
            m_select_dof * vertex_matrix_to_dof_matrix(displacement_map, dim);
        m_displacement_dof_map.makeCompressed();
    }

    // ========================================================================

    // Set vertices at rest using full → reduced map
    m_rest_positions = m_select_vertices * full_rest_positions;
    // m_rest_positions = vertices(full_rest_positions);

    // Map faces and edges to only included vertices
    if (!include_all_vertices) {
        for (int i = 0; i < m_edges.rows(); i++) {
            for (int j = 0; j < m_edges.cols(); j++) {
                long new_id = m_full_vertex_to_vertex[m_edges(i, j)];
                assert(new_id >= 0 && new_id < num_vertices());
                m_edges(i, j) = new_id;
            }
        }

        for (int i = 0; i < m_faces.rows(); i++) {
            for (int j = 0; j < m_faces.cols(); j++) {
                long new_id = m_full_vertex_to_vertex[m_faces(i, j)];
                assert(new_id >= 0 && new_id < num_vertices());
                m_faces(i, j) = new_id;
            }
        }
    } // else no need to change the edges and faces

    m_faces_to_edges = construct_faces_to_edges(m_faces, m_edges);
    construct_edges_to_faces();

    init_codim_vertices();
    init_codim_edges();
    init_areas();
    init_adjacencies();
    // Compute these manually if needed.
    // init_area_jacobian();
}

// ============================================================================

void CollisionMesh::construct_edges_to_faces()
{
    if (dim() == 2)
        return;
    
    m_edges_to_faces.setOnes(num_edges(), 2);
    m_edges_to_faces *= -1;
    for (int f = 0; f < m_faces_to_edges.rows(); f++)
    {
        for (int le = 0; le < 3; le++)
        {
            if (m_edges_to_faces(m_faces_to_edges(f, le), 0) < 0)
                m_edges_to_faces(m_faces_to_edges(f, le), 0) = f;
            else if (m_edges_to_faces(m_faces_to_edges(f, le), 1) < 0)
                m_edges_to_faces(m_faces_to_edges(f, le), 1) = f;
            else
                assert(false);
        }
    }
}

void CollisionMesh::init_codim_vertices()
{
    m_is_codim_vertex.assign(num_vertices(), true);
    for (int i : m_edges.reshaped()) {
        m_is_codim_vertex[i] = false;
    }

    m_codim_vertices.resize(
        std::count(m_is_codim_vertex.begin(), m_is_codim_vertex.end(), true));

    int j = 0;
    for (int i = 0; i < num_vertices(); i++) {
        if (m_is_codim_vertex[i]) {
            assert(j < m_codim_vertices.size());
            m_codim_vertices[j++] = i;
        }
    }
    assert(j == m_codim_vertices.size());
}

void CollisionMesh::init_codim_edges()
{
    m_is_codim_edge.assign(num_edges(), true);
    for (int i : m_faces_to_edges.reshaped()) {
        m_is_codim_edge[i] = false;
    }

    m_codim_edges.resize(
        std::count(m_is_codim_edge.begin(), m_is_codim_edge.end(), true));

    int j = 0;
    for (int i = 0; i < num_edges(); i++) {
        if (m_is_codim_edge[i]) {
            assert(j < m_codim_edges.size());
            m_codim_edges[j++] = i;
        }
    }
    assert(j == m_codim_edges.size());
}

// ============================================================================

void CollisionMesh::init_selection_matrices(const int dim)
{
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(num_vertices());
    for (int vi = 0; vi < num_vertices(); vi++) {
        triplets.emplace_back(vi, m_vertex_to_full_vertex[vi], 1.0);
    }

    m_select_vertices.resize(num_vertices(), full_num_vertices());
    m_select_vertices.setFromTriplets(triplets.begin(), triplets.end());
    m_select_vertices.makeCompressed();

    m_select_dof = vertex_matrix_to_dof_matrix(m_select_vertices, dim);
}

Eigen::SparseMatrix<double> CollisionMesh::vertex_matrix_to_dof_matrix(
    const Eigen::SparseMatrix<double>& M_V, int dim)
{
    std::vector<Eigen::Triplet<double>> triplets;
    using InnerIterator = Eigen::SparseMatrix<double>::InnerIterator;
    for (int k = 0; k < M_V.outerSize(); ++k) {
        for (InnerIterator it(M_V, k); it; ++it) {
            for (int d = 0; d < dim; d++) {
                triplets.emplace_back(
                    dim * it.row() + d, dim * it.col() + d, it.value());
            }
        }
    }

    Eigen::SparseMatrix<double> M_dof(M_V.rows() * dim, M_V.cols() * dim);
    M_dof.setFromTriplets(triplets.begin(), triplets.end());
    M_dof.makeCompressed();
    return M_dof;
}

// ============================================================================

void CollisionMesh::init_adjacencies()
{
    m_vertex_vertex_adjacencies.resize(num_vertices());
    m_vertex_edge_adjacencies.resize(num_vertices());
    // Edges includes the edges of the faces
    for (int i = 0; i < m_edges.rows(); i++) {
        m_vertex_vertex_adjacencies[m_edges(i, 0)].insert(m_edges(i, 1));
        m_vertex_vertex_adjacencies[m_edges(i, 1)].insert(m_edges(i, 0));
        m_vertex_edge_adjacencies[m_edges(i, 0)].insert(i);
        m_vertex_edge_adjacencies[m_edges(i, 1)].insert(i);
    }

    m_edge_vertex_adjacencies.resize(m_edges.rows());
    for (int i = 0; i < m_faces.rows(); i++) {
        for (int j = 0; j < 3; ++j) {
            m_edge_vertex_adjacencies[m_faces_to_edges(i, j)].insert(
                m_faces(i, (j + 2) % 3));
        }
    }

    // Is the vertex on the boundary of the triangle mesh in 3D or polyline in
    // 2D
    m_is_vertex_on_boundary.resize(num_vertices(), true);
    if (dim() == 2) {
        for (int i = 0; i < num_vertices(); i++) {
            m_is_vertex_on_boundary[i] =
                m_vertex_vertex_adjacencies[i].size() <= 1;
        }
    } else {
        for (int i = 0; i < m_edges.rows(); i++) {
            // If edge is part of two triangles
            if (m_edge_vertex_adjacencies[i].size() >= 2) {
                for (int j = 0; j < 2; j++) {
                    m_is_vertex_on_boundary[m_edges(i, j)] = false;
                }
            }
        }
    }
}

void CollisionMesh::init_areas()
{
    // m_vertices_to_edges.resize(num_vertices());
    // for (int i = 0; i < m_edges.rows(); i++) {
    //     for (int j = 0; j < m_edges.cols(); j++) {
    //         m_vertices_to_edges[m_edges(i, j)].push_back(i);
    //     }
    // }
    //
    m_vertices_to_faces.resize(num_vertices());
    for (int i = 0; i < m_faces.rows(); i++) {
        for (int j = 0; j < m_faces.cols(); j++) {
            m_vertices_to_faces[m_faces(i, j)].push_back(i);
        }
    }

    // Compute vertex areas as the sum of ½ the length of connected edges
    Eigen::VectorXd vertex_edge_areas =
        Eigen::VectorXd::Constant(num_vertices(), -1);
    for (int i = 0; i < m_edges.rows(); i++) {
        const VectorMax3d e0 = m_rest_positions.row(m_edges(i, 0));
        const VectorMax3d e1 = m_rest_positions.row(m_edges(i, 1));
        double edge_len = (e1 - e0).norm();

        for (int j = 0; j < m_edges.cols(); j++) {
            if (vertex_edge_areas[m_edges(i, j)] < 0) {
                vertex_edge_areas[m_edges(i, j)] = 0;
            }
            vertex_edge_areas[m_edges(i, j)] += edge_len / 2;
        }
    }

    // Compute vertex/edge areas as the sum of ⅓ the area of connected face
    Eigen::VectorXd vertex_face_areas =
        Eigen::VectorXd::Constant(num_vertices(), -1);
    m_edge_areas.setConstant(m_edges.rows(), -1);
    if (dim() == 3) {
        for (int i = 0; i < m_faces.rows(); i++) {
            const Eigen::Vector3d f0 = m_rest_positions.row(m_faces(i, 0));
            const Eigen::Vector3d f1 = m_rest_positions.row(m_faces(i, 1));
            const Eigen::Vector3d f2 = m_rest_positions.row(m_faces(i, 2));
            double face_area = (f1 - f0).cross(f2 - f0).norm() / 2;

            for (int j = 0; j < m_faces.cols(); ++j) {
                if (vertex_face_areas[m_faces(i, j)] < 0) {
                    vertex_face_areas[m_faces(i, j)] = 0;
                }
                vertex_face_areas[m_faces(i, j)] += face_area / 3;

                if (m_edge_areas[m_faces_to_edges(i, j)] < 0) {
                    m_edge_areas[m_faces_to_edges(i, j)] = 0;
                }
                m_edge_areas[m_faces_to_edges(i, j)] += face_area / 3;
            }
        }
    }

    // Select the area based on the order face, edge, codim
    m_vertex_areas =
        (vertex_face_areas.array() < 0)
            .select(
                (vertex_edge_areas.array() < 0).select(1, vertex_edge_areas),
                vertex_face_areas);

    for (int i = 0; i < m_edge_areas.size(); i++) {
        if (m_edge_areas[i] < 0) {
            // Use the edge length for codim edges
            const VectorMax3d e0 = m_rest_positions.row(m_edges(i, 0));
            const VectorMax3d e1 = m_rest_positions.row(m_edges(i, 1));
            m_edge_areas[i] = (e1 - e0).norm();
        }
    }
}

void CollisionMesh::init_area_jacobians()
{
    std::vector<bool> was_vertex_visited(num_vertices(), false);
    std::vector<bool> was_edge_visited(num_edges(), false);

    m_vertex_area_jacobian.resize(
        num_vertices(), Eigen::SparseVector<double>(ndof()));
    m_edge_area_jacobian.resize(
        num_edges(), Eigen::SparseVector<double>(ndof()));

    // Compute vertex/edge areas as the sum of ⅓ the area of connected face
    for (int i = 0; i < m_faces.rows(); i++) {
        assert(dim() == 3);
        const Eigen::Vector3d f0 = m_rest_positions.row(m_faces(i, 0));
        const Eigen::Vector3d f1 = m_rest_positions.row(m_faces(i, 1));
        const Eigen::Vector3d f2 = m_rest_positions.row(m_faces(i, 2));

        const Vector9d face_area_gradient =
            triangle_area_gradient(f0, f1, f2) / 3.0;

        for (int j = 0; j < m_faces.cols(); ++j) {
            // compute gradient of area

            was_vertex_visited[m_faces(i, j)] = true;
            local_gradient_to_global_gradient(
                face_area_gradient, m_faces.row(i), dim(),
                m_vertex_area_jacobian[m_faces(i, j)]);

            was_edge_visited[m_faces_to_edges(i, j)] = true;
            local_gradient_to_global_gradient(
                face_area_gradient, m_faces.row(i), dim(),
                m_edge_area_jacobian[m_faces_to_edges(i, j)]);
        }
    }

    // Compute unvisited vertex areas as the sum of ½ the length of connected
    // edges
    for (int i = 0; i < m_edges.rows(); i++) {
        const int e0i = m_edges(i, 0), e1i = m_edges(i, 1);
        const VectorMax3d e0 = m_rest_positions.row(e0i);
        const VectorMax3d e1 = m_rest_positions.row(e1i);

        assert(was_vertex_visited[e0i] == was_vertex_visited[e1i]);
        if (was_vertex_visited[e0i] && was_edge_visited[i]) {
            continue;
        }

        const VectorMax6d edge_len_gradient = edge_length_gradient(e0, e1);

        if (!was_vertex_visited[e0i]) {
            for (int j = 0; j < m_edges.cols(); j++) {
                local_gradient_to_global_gradient(
                    edge_len_gradient / 2, m_edges.row(i), dim(),
                    m_vertex_area_jacobian[m_edges(i, j)]);
            }
        }

        if (!was_edge_visited[i]) {
            local_gradient_to_global_gradient(
                edge_len_gradient, m_edges.row(i), dim(),
                m_edge_area_jacobian[i]);
        }
    }
}

// ============================================================================/

Eigen::MatrixXd
CollisionMesh::vertices(const Eigen::MatrixXd& full_positions) const
{
    // full_U = full_V - full_V_rest
    assert(full_positions.rows() == full_num_vertices());
    assert(full_positions.cols() == dim());
    return displace_vertices(full_positions - m_full_rest_positions);
}

Eigen::MatrixXd CollisionMesh::displace_vertices(
    const Eigen::MatrixXd& full_displacements) const
{
    // V_rest + S * T * full_U; m_displacement_map = S * T
    return m_rest_positions + map_displacements(full_displacements);
}

Eigen::MatrixXd CollisionMesh::map_displacements(
    const Eigen::MatrixXd& full_displacements) const
{
    assert(m_displacement_map.cols() == full_displacements.rows());
    assert(full_displacements.cols() == dim());
    return m_displacement_map * full_displacements;
}

// ============================================================================/

Eigen::VectorXd CollisionMesh::to_full_dof(const Eigen::VectorXd& x) const
{
    // ∇_{full} f(S * T * x_full) = Tᵀ * Sᵀ * ∇_{collision} f(S * T * x_full)
    // x = ∇_{collision} f(S * T * x_full); m_displacement_dof_map = S * T
    return m_displacement_dof_map.transpose() * x;
}

Eigen::SparseMatrix<double>
CollisionMesh::to_full_dof(const Eigen::SparseMatrix<double>& X) const
{
    // ∇_{full} Tᵀ * Sᵀ * ∇_{collision} f(S * T * x_full)
    //      = Tᵀ * Sᵀ * [∇_{collision}² f(S * T * x_full)] * S * T
    // X = ∇_{collision}² f(S * T * x_full); m_displacement_dof_map = S * T
    return m_displacement_dof_map.transpose() * X * m_displacement_dof_map;
}

// ============================================================================/

std::vector<bool> CollisionMesh::construct_is_on_surface(
    const long num_vertices,
    const Eigen::MatrixXi& edges,
    const Eigen::VectorXi& codim_vertices)
{
    std::vector<bool> is_on_surface(num_vertices, false);
    for (int i = 0; i < codim_vertices.size(); i++) {
        assert(codim_vertices[i] < num_vertices);
        is_on_surface[codim_vertices[i]] = true;
    }
    // Column first because colmajor
    for (size_t ej = 0; ej < edges.cols(); ej++) {
        for (size_t ei = 0; ei < edges.rows(); ei++) {
            assert(edges(ei, ej) < num_vertices);
            is_on_surface[edges(ei, ej)] = true;
        }
    }
    return is_on_surface;
}

// ============================================================================/

Eigen::MatrixXi CollisionMesh::construct_faces_to_edges(
    const Eigen::MatrixXi& faces, const Eigen::MatrixXi& edges)
{
    if (faces.size() == 0) {
        return Eigen::MatrixXi(faces.rows(), faces.cols());
    }
    assert(edges.size() != 0);

    unordered_map<std::pair<int, int>, size_t> edge_map;
    for (size_t ei = 0; ei < edges.rows(); ei++) {
        edge_map.emplace(
            std::make_pair<int, int>(
                edges.row(ei).minCoeff(), edges.row(ei).maxCoeff()),
            ei);
    }

    Eigen::MatrixXi faces_to_edges(faces.rows(), faces.cols());
    for (int fi = 0; fi < faces.rows(); fi++) {
        for (int fj = 0; fj < faces.cols(); fj++) {
            const int vi = faces(fi, fj);
            const int vj = faces(fi, (fj + 1) % faces.cols());
            std::pair<int, int> e(std::min(vi, vj), std::max(vi, vj));
            auto search = edge_map.find(e);
            if (search != edge_map.end()) {
                faces_to_edges(fi, fj) = search->second;
            } else {
                throw std::runtime_error("Unable to find edge!");
            }
        }
    }

    return faces_to_edges;
}

double CollisionMesh::edge_length(const int &edge_id) const
{
    return (m_rest_positions.row(m_edges(edge_id, 0)) - m_rest_positions.row(m_edges(edge_id, 1))).norm();
}

double CollisionMesh::max_edge_length() const
{
    double val = 0;
    for (int i = 0; i < m_edges.rows(); i++)
        val = std::max(edge_length(i), val);
    return val;
}

std::vector<long> CollisionMesh::find_vertex_adjacent_vertices(const long &v) const
{
    std::vector<long> neighbors;
    if (dim() == 2)
    {
        neighbors.assign(2, -1);
        for (long i : vertex_edge_adjacencies()[v])
        {
            if (edges()(i, 0) == v)
                neighbors[0] = edges()(i, 1);
            else if (edges()(i, 1) == v)
                neighbors[1] = edges()(i, 0);
            else
                throw std::runtime_error("Invalid edge-vertex adjacency!");
        }
    }
    else
    {
        std::unordered_map<long, long> map;
        for (auto f : vertices_to_faces()[v])
        {
            for (int lv = 0; lv < 3; lv++)
            {
                if (faces()(f, lv) == v)
                {
                    map[faces()(f, (lv+1)%3)] = faces()(f, (lv+2)%3);
                    break;
                }
            }
        }
        if (vertices_to_faces()[v].size() != map.size())
            throw std::runtime_error("Non-manifold vertex! Map size smaller than neighbor!");
        
        auto iter = map.find(map.begin()->first);
        while (neighbors.empty() || iter->first != neighbors.front())
        {
            neighbors.push_back(iter->first);
            iter = map.find(iter->second);
            if (iter == map.end())
            {
                logger().error("neighbor faces {}, map {}", vertices_to_faces()[v].size(), map);
                throw std::runtime_error("Non-manifold vertex! Cannot find next neighbor!");
            }
        }
        if (neighbors.size() != map.size())
            throw std::runtime_error("Non-manifold vertex!");
    }
    return neighbors;
}

/// @brief 
/// @param e edge id 
/// @return 4 vertices of the 2 faces adjacent to this edge, in the order of [e0, e1, f0, f1]
/// the two faces are [f0, e0, e1] and [f1, e1, e0]
std::array<long, 4> CollisionMesh::find_edge_adjacent_vertices(const long &e) const
{
    std::array<long, 4> neighbors;
    if (dim() == 2)
        throw std::runtime_error("find_edge_adjacent_vertices is only for 3d collision mesh!");
    
    std::array<long, 2> lf = {{edges_to_faces()(e, 0), edges_to_faces()(e, 1)}};

    for (int j : {0, 1})
    {
        int i;
        for (i = 0; i < 3; i++)
        {
            const auto va = faces()(lf[j], i);

            if (va != edges()(e, 0) && va != edges()(e, 1))
            {
                neighbors[2 + j] = va;
                break;
            }
        }

        if (j == 0)
        {
            neighbors[0] = faces()(lf[j], (i+1) % 3);
            neighbors[1] = faces()(lf[j], (i+2) % 3);
        }
    }

    return neighbors;
}

} // namespace ipc
