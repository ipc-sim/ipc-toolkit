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

    size_t num_vertices() const { return vertices_at_rest().rows(); }
    size_t dim() const { return vertices_at_rest().cols(); }
    size_t ndof() const { return num_vertices() * dim(); }

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
    void init_dof_to_full_dof();

    size_t full_ndof() const { return full_vertex_to_vertex.size() * dim(); }

    Eigen::MatrixXd m_vertices_at_rest;
    /// Edges as rows of indicies into V.
    Eigen::MatrixXi m_edges;
    /// Triangular faces as rows of indicies into V.
    Eigen::MatrixXi m_faces;
    /// Map from F edges to rows of E.
    Eigen::MatrixXi m_faces_to_edges;

    Eigen::VectorXi full_vertex_to_vertex;
    Eigen::VectorXi vertex_to_full_vertex;
    Eigen::VectorXi dof_to_full_dof;
};

} // namespace ipc
