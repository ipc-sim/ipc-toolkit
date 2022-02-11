#include <ipc/collision_mesh.hpp>

#include <ipc/utils/unordered_map_and_set.hpp>

#include <igl/slice.h>
#include <igl/slice_into.h>

namespace ipc {

CollisionMesh::CollisionMesh(
    const Eigen::MatrixXd& surface_vertices_at_rest,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces)
    : m_vertices_at_rest(surface_vertices_at_rest)
    , m_edges(edges)
    , m_faces(faces)
{
    m_faces_to_edges = construct_faces_to_edges(m_faces, m_edges);

    // Assumes surface is full mesh
    m_full_size = surface_vertices_at_rest.rows();
    m_full_to_surface.setLinSpaced(full_size(), 0, full_size() - 1);
    assert(full_to_surface().size() == full_size());
    m_surface_to_full = m_full_to_surface;
}

CollisionMesh::CollisionMesh(
    const std::vector<bool>& is_on_surface,
    const Eigen::MatrixXd& full_vertices_at_rest,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces)
    : m_edges(edges)
    , m_faces(faces)
{
    // Assumes surface is full mesh
    m_full_size = full_vertices_at_rest.rows();
    assert(is_on_surface.size() == full_size());

    m_full_to_surface.setConstant(full_size(), -1);
    std::vector<int> dynamic_surface_to_full;
    size_t surface_size = 0;
    for (size_t i = 0; i < full_size(); i++) {
        if (is_on_surface[i]) {
            m_full_to_surface(i) = surface_size;
            dynamic_surface_to_full.push_back(i);
            surface_size++;
        }
    }
    m_surface_to_full = Eigen::Map<Eigen::VectorXi>(
        dynamic_surface_to_full.data(), dynamic_surface_to_full.size());

    m_vertices_at_rest = surface_vertices(full_vertices_at_rest);

    for (int i = 0; i < m_edges.rows(); i++) {
        for (int j = 0; j < m_edges.cols(); j++) {
            long new_id = full_to_surface()[m_edges(i, j)];
            assert(new_id >= 0 && new_id < surface_size);
            m_edges(i, j) = new_id;
        }
    }

    for (int i = 0; i < m_faces.rows(); i++) {
        for (int j = 0; j < m_faces.cols(); j++) {
            long new_id = full_to_surface()[m_faces(i, j)];
            assert(new_id >= 0 && new_id < surface_size);
            m_faces(i, j) = new_id;
        }
    }

    m_faces_to_edges = construct_faces_to_edges(m_faces, m_edges);
}

////////////////////////////////////////////////////////////////////////////////

Eigen::VectorXd
CollisionMesh::map_surface_to_full(const Eigen::VectorXd& x) const
{
    Eigen::VectorXd full_x = Eigen::VectorXd::Zero(full_size() * dim());
    Eigen::VectorXi surface_to_full_flat(x.size());
    for (int i = 0; i < surface_size(); i++) {
        for (int d = 0; d < dim(); d++) {
            surface_to_full_flat[dim() * i + d] =
                dim() * surface_to_full()[i] + d;
        }
    }
    igl::slice_into(x, surface_to_full_flat, full_x);
    return full_x;
}

Eigen::SparseMatrix<double>
CollisionMesh::map_surface_to_full(const Eigen::SparseMatrix<double>& X) const
{
    // initializes to zero
    int n = full_size() * dim();
    Eigen::SparseMatrix<double> full_X(n, n);
    full_X.reserve(X.nonZeros());

    Eigen::VectorXi surface_to_full_flat(surface_size() * dim());
    for (int i = 0; i < surface_size(); i++) {
        for (int d = 0; d < dim(); d++) {
            surface_to_full_flat[dim() * i + d] =
                dim() * surface_to_full()[i] + d;
        }
    }

    igl::slice_into(X, surface_to_full_flat, surface_to_full_flat, full_X);
    full_X.makeCompressed();
    return full_X;
}

Eigen::MatrixXd
CollisionMesh::surface_vertices(const Eigen::MatrixXd& full_V) const
{
    Eigen::MatrixXd surface_V(surface_size(), full_V.cols());
    igl::slice(full_V, surface_to_full(), /*dim=*/1, surface_V);
    return surface_V;
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