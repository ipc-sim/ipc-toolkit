#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>

namespace ipc {

class CollisionMesh {
public:
    CollisionMesh() { }

    CollisionMesh(
        const Eigen::MatrixXd& vertices_at_rest,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces);

    CollisionMesh(
        const std::vector<bool>& include_vertex,
        const Eigen::MatrixXd& full_vertices_at_rest,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces);

    ~CollisionMesh() { }

    const Eigen::MatrixXd& vertices_at_rest() const
    {
        return m_vertices_at_rest;
    }
    Eigen::MatrixXd vertices(const Eigen::MatrixXd& full_V) const;
    const Eigen::MatrixXi& edges() const { return m_edges; }
    const Eigen::MatrixXi& faces() const { return m_faces; }
    const Eigen::MatrixXi& faces_to_edges() const { return m_faces_to_edges; }

    void
    set_linear_vertex_map(const Eigen::SparseMatrix<double>& linear_vertex_map);

    size_t num_vertices() const { return m_num_vertices; }
    size_t dim() const { return m_dim; }
    size_t ndof() const { return num_vertices() * dim(); }
    size_t full_num_vertices() const { return m_full_num_vertices; }
    size_t full_ndof() const { return full_num_vertices() * dim(); }

    size_t to_full_vertex_id(const size_t id)
    {
        assert(id < num_vertices());
        return vertex_to_full_vertex[id];
    }
    Eigen::VectorXd to_full_dof(const Eigen::VectorXd& x) const;
    Eigen::SparseMatrix<double>
    to_full_dof(const Eigen::SparseMatrix<double>& X) const;

    static std::vector<bool> construct_is_on_surface(
        const int num_vertices, const Eigen::MatrixXi& edges);

    static CollisionMesh build_from_full_mesh(
        const Eigen::MatrixXd& full_vertices_at_rest,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces)
    {
        return CollisionMesh(
            construct_is_on_surface(full_vertices_at_rest.rows(), edges),
            full_vertices_at_rest, edges, faces);
    }

    /// Construct a matrix that maps from the faces' edges to rows in the edges
    /// matrix.
    static Eigen::MatrixXi construct_faces_to_edges(
        const Eigen::MatrixXi& faces, const Eigen::MatrixXi& edges);

    /// A function that takes two vertex IDs (row numbers in V) and returns true
    /// if the vertices (and faces or edges containing the vertices) can
    /// collide. By default all primitives can collide with all other
    /// primitives.
    std::function<bool(size_t, size_t)> can_collide = [](size_t, size_t) {
        return true;
    };

protected:
    void init_selection_matrix();
    void set_identity_linear_vertex_map();

    // Convert a matrix meant for M_V * V to M_dof * x by duplicating the
    // entries dim times.
    static Eigen::SparseMatrix<double> vertex_matrix_to_dof_matrix(
        const Eigen::SparseMatrix<double>& M_V, int dim);

    Eigen::MatrixXd m_vertices_at_rest;
    /// Edges as rows of indicies into V.
    Eigen::MatrixXi m_edges;
    /// Triangular faces as rows of indicies into V.
    Eigen::MatrixXi m_faces;
    /// Map from F edges to rows of E.
    Eigen::MatrixXi m_faces_to_edges;

    Eigen::VectorXi full_vertex_to_vertex;
    Eigen::VectorXi vertex_to_full_vertex;

    // Selection matrix S ∈ ℝ^{collision×full}
    Eigen::SparseMatrix<double> select_vertices;
    Eigen::SparseMatrix<double> select_dof;

    /// Mapping from full vertices to collision vertices
    Eigen::SparseMatrix<double> m_full_to_collision_vertices;
    /// Mapping from collision DOF to full DOF
    Eigen::SparseMatrix<double> m_full_to_collision_dof;

    int m_full_num_vertices;
    int m_num_vertices;
    int m_dim;
};

} // namespace ipc
