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
    m_faces_to_edges = construct_faces_to_edges(m_faces, m_edges);

    // Assumes collision mesh is full mesh
    full_vertex_to_vertex.setLinSpaced(num_vertices(), 0, num_vertices() - 1);
    vertex_to_full_vertex = full_vertex_to_vertex;
    init_dof_to_full_dof();
}

CollisionMesh::CollisionMesh(
    const std::vector<bool>& include_vertex,
    const Eigen::MatrixXd& full_vertices_at_rest,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces)
    : m_edges(edges)
    , m_faces(faces)
{
    size_t num_full_verts = full_vertices_at_rest.rows();
    assert(include_vertex.size() == num_full_verts);

    full_vertex_to_vertex.setConstant(num_full_verts, -1);
    std::vector<int> dynamic_vertex_to_full_vertex;
    size_t num_verts = 0;
    for (size_t i = 0; i < num_full_verts; i++) {
        if (include_vertex[i]) {
            full_vertex_to_vertex[i] = num_verts;
            dynamic_vertex_to_full_vertex.push_back(i);
            num_verts++;
        }
    }
    vertex_to_full_vertex = Eigen::Map<Eigen::VectorXi>(
        dynamic_vertex_to_full_vertex.data(),
        dynamic_vertex_to_full_vertex.size());

    m_vertices_at_rest = vertices(full_vertices_at_rest);

    for (int i = 0; i < m_edges.rows(); i++) {
        for (int j = 0; j < m_edges.cols(); j++) {
            long new_id = full_vertex_to_vertex[m_edges(i, j)];
            assert(new_id >= 0 && new_id < num_verts);
            m_edges(i, j) = new_id;
        }
    }

    for (int i = 0; i < m_faces.rows(); i++) {
        for (int j = 0; j < m_faces.cols(); j++) {
            long new_id = full_vertex_to_vertex[m_faces(i, j)];
            assert(new_id >= 0 && new_id < num_verts);
            m_faces(i, j) = new_id;
        }
    }

    m_faces_to_edges = construct_faces_to_edges(m_faces, m_edges);

    init_dof_to_full_dof();
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
    Eigen::VectorXd full_x = Eigen::VectorXd::Zero(full_ndof());
    igl::slice_into(x, dof_to_full_dof, full_x);
    return full_x;
}

Eigen::SparseMatrix<double>
CollisionMesh::to_full_dof(const Eigen::SparseMatrix<double>& X) const
{
    // initializes to zero
    Eigen::SparseMatrix<double> full_X(full_ndof(), full_ndof());
    full_X.reserve(X.nonZeros());
    igl::slice_into(X, dof_to_full_dof, dof_to_full_dof, full_X);
    full_X.makeCompressed();
    return full_X;
}

Eigen::MatrixXd CollisionMesh::vertices(const Eigen::MatrixXd& full_V) const
{
    Eigen::MatrixXd V(num_vertices(), full_V.cols());
    igl::slice(full_V, vertex_to_full_vertex, /*dim=*/1, V);
    return V;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<bool> CollisionMesh::construct_is_on_surface(
    const int num_vertices, const Eigen::MatrixXi& edges)
{
    std::vector<bool> include_vertex(num_vertices, false);
    // Column first because colmajor
    for (size_t ej = 0; ej < edges.cols(); ej++) {
        for (size_t ei = 0; ei < edges.rows(); ei++) {
            assert(edges(ei, ej) < num_vertices);
            include_vertex[edges(ei, ej)] = true;
        }
    }
    return include_vertex;
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