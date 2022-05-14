#include <ipc/collision_mesh.hpp>

#include <ipc/utils/unordered_map_and_set.hpp>
#include <ipc/utils/logger.hpp>

namespace ipc {

CollisionMesh::CollisionMesh(
    const Eigen::MatrixXd& vertices_at_rest,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces)
    : CollisionMesh(
        std::vector<bool>(vertices_at_rest.rows(), true),
        vertices_at_rest,
        edges,
        faces)
{
}

bool all(const std::vector<bool>& v)
{
    for (bool b : v) {
        if (!b) {
            return false;
        }
    }
    return true;
}

CollisionMesh::CollisionMesh(
    const std::vector<bool>& include_vertex,
    const Eigen::MatrixXd& full_vertices_at_rest,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces)
    : m_edges(edges)
    , m_faces(faces)
{
    assert(include_vertex.size() == full_vertices_at_rest.rows());
    bool include_all_vertices = all(include_vertex);

    m_full_num_vertices = full_vertices_at_rest.rows();
    m_dim = full_vertices_at_rest.cols();

    if (include_all_vertices) {
        m_num_vertices = full_num_vertices();
        // set full ↔ reduced ≡ identity
        full_vertex_to_vertex.setLinSpaced(
            num_vertices(), 0, num_vertices() - 1);
        vertex_to_full_vertex = full_vertex_to_vertex;
    } else {
        full_vertex_to_vertex.setConstant(full_num_vertices(), -1);
        std::vector<int> dynamic_vertex_to_full_vertex;
        m_num_vertices = 0;
        for (size_t i = 0; i < full_num_vertices(); i++) {
            if (include_vertex[i]) {
                full_vertex_to_vertex[i] = m_num_vertices;
                dynamic_vertex_to_full_vertex.push_back(i);
                m_num_vertices++;
            }
        }
        vertex_to_full_vertex = Eigen::Map<Eigen::VectorXi>(
            dynamic_vertex_to_full_vertex.data(),
            dynamic_vertex_to_full_vertex.size());
    }
    init_selection_matrix();
    set_identity_linear_vertex_map();

    // Set vertices at rest using full → reduced map
    m_vertices_at_rest = m_select_vertices * m_full_vertices_at_rest;
    // m_vertices_at_rest = vertices(full_vertices_at_rest);

    // Map faces and edges to only included vertices
    if (!include_all_vertices) {
        for (int i = 0; i < m_edges.rows(); i++) {
            for (int j = 0; j < m_edges.cols(); j++) {
                long new_id = full_vertex_to_vertex[m_edges(i, j)];
                assert(new_id >= 0 && new_id < m_num_vertices);
                m_edges(i, j) = new_id;
            }
        }

        for (int i = 0; i < m_faces.rows(); i++) {
            for (int j = 0; j < m_faces.cols(); j++) {
                long new_id = full_vertex_to_vertex[m_faces(i, j)];
                assert(new_id >= 0 && new_id < m_num_vertices);
                m_faces(i, j) = new_id;
            }
        }
    } // else no need to change the edges and faces

    m_faces_to_edges = construct_faces_to_edges(m_faces, m_edges);
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
    return M_dof;
}

void CollisionMesh::init_selection_matrix()
{
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(num_vertices());
    for (int vi = 0; vi < num_vertices(); vi++) {
        triplets.emplace_back(vi, vertex_to_full_vertex[vi], 1.0);
    }

    m_select_vertices.resize(num_vertices(), full_num_vertices());
    m_select_vertices.setFromTriplets(triplets.begin(), triplets.end());
    m_select_vertices.makeCompressed();

    m_select_dof = vertex_matrix_to_dof_matrix(m_select_vertices, dim());
    m_select_dof.makeCompressed();
}

void CollisionMesh::set_identity_linear_vertex_map()
{
    // Initilize linear map with identity
    // S * T = S * I = S
    m_full_to_collision_vertices = m_select_vertices;
    m_full_to_collision_dof = m_select_dof;
}

void CollisionMesh::set_linear_vertex_map(
    const Eigen::SparseMatrix<double>& T_vertices)
{
    assert(T_vertices.rows() == full_num_vertices());

    m_full_to_collision_vertices = m_select_vertices * T_vertices;
    m_full_to_collision_vertices.makeCompressed();

    m_full_to_collision_dof =
        m_select_dof * vertex_matrix_to_dof_matrix(T_vertices, dim());
    m_full_to_collision_dof.makeCompressed();
}

////////////////////////////////////////////////////////////////////////////////

// Eigen::MatrixXd CollisionMesh::vertices(const Eigen::MatrixXd& full_V) const
// {
//     // full_U = full_V - full_V_rest
//     // assert(full_V.rows() == full_num_vertices());
//     // assert(full_V.cols() == dim());
//     // return vertices_from_displacements(full_V - m_full_vertices_at_rest);
//     // S * full_V; m_select_vertices = S
//     assert(full_V.rows() == m_select_vertices.cols());
//     assert(full_V.cols() == dim());
//     return m_select_vertices * full_V;
// }

Eigen::MatrixXd
CollisionMesh::vertices_from_displacements(const Eigen::MatrixXd& full_U) const
{
    // V_rest + S * T * full_U; m_full_to_collision_vertices = S * T
    assert(full_U.rows() == m_full_to_collision_vertices.cols());
    assert(full_U.cols() == dim());
    return m_vertices_at_rest + m_full_to_collision_vertices * full_U;
}

Eigen::VectorXd CollisionMesh::to_full_dof(const Eigen::VectorXd& x) const
{
    // ∇_{full} f(S * T * x_full) = Tᵀ * Sᵀ * ∇_{collision} f(S * T * x_full)
    // x = ∇_{collision} f(S * T * x_full); m_full_to_collision_dof = S * T
    return m_full_to_collision_dof.transpose() * x;
}

Eigen::SparseMatrix<double>
CollisionMesh::to_full_dof(const Eigen::SparseMatrix<double>& X) const
{
    // ∇_{full} Tᵀ * Sᵀ * ∇_{collision} f(S * T * x_full)
    //      = Tᵀ * Sᵀ * ∇_{collision}² f(S * T * x_full) * S * T
    // X = ∇_{collision}² f(S * T * x_full); m_full_to_collision_dof = S * T
    return m_full_to_collision_dof.transpose() * X * m_full_to_collision_dof;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<bool> CollisionMesh::construct_is_on_surface(
    const int num_vertices, const Eigen::MatrixXi& edges)
{
    std::vector<bool> is_on_surface(num_vertices, false);
    // Column first because colmajor
    for (size_t ej = 0; ej < edges.cols(); ej++) {
        for (size_t ei = 0; ei < edges.rows(); ei++) {
            assert(edges(ei, ej) < num_vertices);
            is_on_surface[edges(ei, ej)] = true;
        }
    }
    return is_on_surface;
}

////////////////////////////////////////////////////////////////////////////////

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

} // namespace ipc
