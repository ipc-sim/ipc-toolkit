#include <ipc/collision_mesh.hpp>

#include <igl/slice.h>
#include <igl/slice_into.h>

#include <ipc/utils/eigen_ext.hpp>

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

    init_adjacencies();
    init_areas();
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

    init_adjacencies();
    init_areas();
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

void CollisionMesh::init_adjacencies()
{
    m_point_point_adjacencies.resize(num_vertices());
    // Edges includes the edges of the faces
    for (int i = 0; i < m_edges.rows(); i++) {
        m_point_point_adjacencies[m_edges(i, 0)].insert(m_edges(i, 1));
        m_point_point_adjacencies[m_edges(i, 1)].insert(m_edges(i, 0));
    }

    m_edge_point_adjacencies.resize(m_edges.rows());
    for (int i = 0; i < m_faces.rows(); i++) {
        for (int j = 0; j < 3; ++j) {
            m_edge_point_adjacencies[m_faces_to_edges(i, j)].insert(
                m_faces(i, (j + 2) % 3));
        }
    }

    // Is the point on the boundary of the triangle mesh in 3D or polyline in 2D
    m_is_point_on_boundary.resize(num_vertices(), true);
    if (dim() == 2) {
        for (int i = 0; i < num_vertices(); i++) {
            m_is_point_on_boundary[i] =
                m_point_point_adjacencies[i].size() <= 1;
        }
    } else {
        for (int i = 0; i < m_edges.rows(); i++) {
            // If edge is part of two triangles
            if (m_edge_point_adjacencies[i].size() >= 2) {
                for (int j = 0; j < 2; j++) {
                    m_is_point_on_boundary[m_edges(i, j)] = false;
                }
            }
        }
    }
}

void CollisionMesh::init_areas()
{
    // Compute point areas as the sum of ½ the length of connected edges
    Eigen::VectorXd point_edge_areas =
        Eigen::VectorXd::Constant(num_vertices(), -1);
    for (int i = 0; i < m_edges.rows(); i++) {
        const auto& e0 = m_vertices_at_rest.row(m_edges(i, 0));
        const auto& e1 = m_vertices_at_rest.row(m_edges(i, 1));
        double edge_len = (e1 - e0).norm();
        for (int j = 0; j < m_edges.cols(); j++) {
            if (point_edge_areas[m_edges(i, j)] < 0) {
                point_edge_areas[m_edges(i, j)] = 0;
            }
            point_edge_areas[m_edges(i, j)] += edge_len / 2;
        }
    }

    // Compute point/edge areas as the sum of ⅓ the area of connected face
    Eigen::VectorXd point_face_areas =
        Eigen::VectorXd::Constant(num_vertices(), -1);
    m_edge_areas.setConstant(m_edges.rows(), -1);
    if (dim() == 3) {
        for (int i = 0; i < m_faces.rows(); i++) {
            const auto& f0 = m_vertices_at_rest.row(m_faces(i, 0));
            const auto& f1 = m_vertices_at_rest.row(m_faces(i, 1));
            const auto& f2 = m_vertices_at_rest.row(m_faces(i, 2));
            double face_area = cross(f1 - f0, f2 - f0).norm() / 2;

            for (int j = 0; j < m_faces.cols(); ++j) {
                if (point_face_areas[m_faces(i, j)] < 0) {
                    point_face_areas[m_faces(i, j)] = 0;
                }
                point_face_areas[m_faces(i, j)] += face_area / 3;

                if (m_edge_areas[m_faces_to_edges(i, j)] < 0) {
                    m_edge_areas[m_faces_to_edges(i, j)] = 0;
                }
                m_edge_areas[m_faces_to_edges(i, j)] += face_area / 3;
            }
        }
    }

    // Select the area based on the order face, edge, codim
    m_point_areas =
        (point_face_areas.array() < 0)
            .select(
                (point_edge_areas.array() < 0).select(1, point_edge_areas),
                point_face_areas);

    // Select the area based on the order face, codim
    m_edge_areas = (m_edge_areas.array() < 0).select(1, m_edge_areas);
}

////////////////////////////////////////////////////////////////////////////////

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
