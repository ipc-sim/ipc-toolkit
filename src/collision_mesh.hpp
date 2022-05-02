#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <ipc/utils/unordered_map_and_set.hpp>

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

    const std::vector<unordered_set<int>>& point_point_adjacencies() const
    {
        return m_point_point_adjacencies;
    }
    const std::vector<unordered_set<int>>& edge_point_adjacencies() const
    {
        return m_edge_point_adjacencies;
    }
    bool is_point_on_boundary(const int i) const
    {
        return m_is_point_on_boundary[i];
    }
    const Eigen::VectorXd& point_areas() const { return m_point_areas; }
    const Eigen::VectorXd& edge_areas() const { return m_edge_areas; }

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
    void init_adjacencies();
    void init_areas();

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

    /// Points adjacent to points
    std::vector<unordered_set<int>> m_point_point_adjacencies;
    /// Edges adjacent to edges
    std::vector<unordered_set<int>> m_edge_point_adjacencies;

    /// Is point on the boundary of the triangle mesh in 3D or polyline in 2D?
    std::vector<bool> m_is_point_on_boundary;

    Eigen::VectorXd m_point_areas;
    Eigen::VectorXd m_edge_areas;
};

} // namespace ipc
