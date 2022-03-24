#include <ipc/collision_mesh.hpp>

#include <ipc/utils/unordered_map_and_set.hpp>

#include <igl/slice.h>
#include <igl/slice_into.h>

namespace ipc {

CollisionMesh::CollisionMesh(
    const Eigen::MatrixXd& vertices_at_rest,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces)
    : m_vertices_at_rest(vertices_at_rest)
    , m_edges(edges)
    , m_faces(faces)
{
    m_num_vertices = m_full_num_vertices = vertices_at_rest.rows();
    m_dim = vertices_at_rest.cols();

    m_faces_to_edges = construct_faces_to_edges(m_faces, m_edges);

    // Assumes collision mesh is full mesh
    full_vertex_to_vertex.setLinSpaced(num_vertices(), 0, num_vertices() - 1);
    vertex_to_full_vertex = full_vertex_to_vertex;
    init_dof_to_full_dof();

    set_identity_linear_vertex_map();
}

CollisionMesh::CollisionMesh(
    const std::vector<bool>& include_vertex,
    const Eigen::MatrixXd& full_vertices_at_rest,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces)
    : m_edges(edges)
    , m_faces(faces)
{
    m_full_num_vertices = full_vertices_at_rest.rows();
    m_dim = full_vertices_at_rest.cols();

    assert(include_vertex.size() == full_num_vertices());

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

    set_identity_linear_vertex_map();

    m_vertices_at_rest = vertices(full_vertices_at_rest);

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

    m_faces_to_edges = construct_faces_to_edges(m_faces, m_edges);

    init_dof_to_full_dof();
}

void CollisionMesh::set_identity_linear_vertex_map()
{
    // Initilize linear map with identity
    m_linear_vertex_map.resize(full_num_vertices(), full_num_vertices());
    m_linear_vertex_map.setIdentity();
    m_linear_dof_map.resize(full_ndof(), full_ndof());
    m_linear_dof_map.setIdentity();
}

void CollisionMesh::set_linear_vertex_map(
    const Eigen::SparseMatrix<double>& linear_vertex_map)
{
    int n = linear_vertex_map.rows(), m = linear_vertex_map.cols();

    assert(n == num_vertices());
    m_linear_vertex_map = linear_vertex_map;

    std::vector<Eigen::Triplet<double>> triplets;
    using InnerIterator = Eigen::SparseMatrix<double>::InnerIterator;
    for (int k = 0; k < m_linear_vertex_map.outerSize(); ++k) {
        for (InnerIterator it(m_linear_vertex_map, k); it; ++it) {
            for (int d = 0; d < dim(); d++) {
                triplets.emplace_back(
                    dim() * it.row() + d, dim() * it.col() + d, it.value());
            }
        }
    }

    m_linear_dof_map.resize(n * dim(), m * dim());
    m_linear_dof_map.setFromTriplets(triplets.begin(), triplets.end());
    m_linear_dof_map.makeCompressed();
}

////////////////////////////////////////////////////////////////////////////////

void CollisionMesh::init_dof_to_full_dof()
{
    dof_to_full_dof.resize(ndof());
    for (int i = 0; i < num_vertices(); i++) {
        for (int d = 0; d < dim(); d++) {
            dof_to_full_dof[dim() * i + d] =
                dim() * vertex_to_full_vertex[i] + d;
        }
    }
}

Eigen::VectorXd CollisionMesh::to_full_dof(const Eigen::VectorXd& x) const
{
    // ∇_{full} f(S * T * x_full) = Tᵀ * Sᵀ * ∇_{collision} f(S * T * x_full)
    // x = ∇_{collision} f(S * T * x_full)
    Eigen::VectorXd full_x = Eigen::VectorXd::Zero(full_ndof());
    igl::slice_into(x, dof_to_full_dof, full_x);
    return m_linear_dof_map.transpose() * full_x;
}

Eigen::SparseMatrix<double>
CollisionMesh::to_full_dof(const Eigen::SparseMatrix<double>& X) const
{
    // initializes to zero
    Eigen::SparseMatrix<double> full_X(full_ndof(), full_ndof());
    full_X.reserve(X.nonZeros());
    igl::slice_into(X, dof_to_full_dof, dof_to_full_dof, full_X);
    full_X.makeCompressed();
    return m_linear_dof_map.transpose() * full_X * m_linear_dof_map;
}

Eigen::MatrixXd CollisionMesh::vertices(const Eigen::MatrixXd& full_V) const
{
    // S * T * full_V
    Eigen::MatrixXd collision_vertices = m_linear_vertex_map * full_V;
    if (vertex_to_full_vertex.size() == collision_vertices.rows())
        return collision_vertices;
    Eigen::MatrixXd V(num_vertices(), dim());
    igl::slice(collision_vertices, vertex_to_full_vertex, /*dim=*/1, V);
    return V;
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

    auto max_vi = edges.maxCoeff();
    auto hash_edge = [&max_vi](const Eigen::RowVector2i& e) {
        return e.minCoeff() * max_vi + e.maxCoeff();
    };
    auto eq_edges = [](const Eigen::RowVector2i& ea,
                       const Eigen::RowVector2i& eb) {
        return ea.minCoeff() == eb.minCoeff() && ea.maxCoeff() == eb.maxCoeff();
    };

    unordered_map<
        Eigen::RowVector2i, int, decltype(hash_edge), decltype(eq_edges)>
        edge_map(/*bucket_count=*/edges.rows(), hash_edge, eq_edges);
    for (int ei = 0; ei < edges.rows(); ei++) {
        edge_map[edges.row(ei)] = ei;
    }

    Eigen::MatrixXi faces_to_edges(faces.rows(), faces.cols());
    for (int fi = 0; fi < faces.rows(); fi++) {
        for (int fj = 0; fj < faces.cols(); fj++) {
            Eigen::RowVector2i e(
                faces(fi, fj), faces(fi, (fj + 1) % faces.cols()));
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
