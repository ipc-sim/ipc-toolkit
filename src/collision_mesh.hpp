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
        const Eigen::MatrixXi& faces,
        const Eigen::SparseMatrix<double>& displacement_map =
            Eigen::SparseMatrix<double>());

    CollisionMesh(
        const std::vector<bool>& include_vertex,
        const Eigen::MatrixXd& full_vertices_at_rest,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        const Eigen::SparseMatrix<double>& displacement_map =
            Eigen::SparseMatrix<double>());

    ~CollisionMesh() { }

    size_t num_vertices() const { return m_vertex_to_full_vertex.size(); }
    size_t num_edges() const { return edges().rows(); }
    size_t num_faces() const { return faces().rows(); }
    size_t dim() const { return m_vertices_at_rest.cols(); }
    size_t ndof() const { return num_vertices() * dim(); }
    size_t full_num_vertices() const { return m_full_vertex_to_vertex.size(); }
    size_t full_ndof() const { return full_num_vertices() * dim(); }

    const Eigen::MatrixXd& vertices_at_rest() const
    {
        return m_vertices_at_rest;
    }
    const Eigen::MatrixXi& edges() const { return m_edges; }
    const Eigen::MatrixXi& faces() const { return m_faces; }
    const Eigen::MatrixXi& faces_to_edges() const { return m_faces_to_edges; }
    // const std::vector<std::vector<int>>& vertices_to_edges() const
    // {
    //     return m_vertices_to_edges;
    // }
    // const std::vector<std::vector<int>>& vertices_to_faces() const
    // {
    //     return m_vertices_to_faces;
    // }

    ///////////////////////////////////////////////////////////////////////////

    Eigen::MatrixXd vertices(const Eigen::MatrixXd& full_vertices) const;

    Eigen::MatrixXd
    displace_vertices(const Eigen::MatrixXd& full_displacements) const;

    size_t to_full_vertex_id(const size_t id) const
    {
        assert(id < num_vertices());
        return m_vertex_to_full_vertex[id];
    }
    Eigen::VectorXd to_full_dof(const Eigen::VectorXd& x) const;
    Eigen::SparseMatrix<double>
    to_full_dof(const Eigen::SparseMatrix<double>& X) const;

    ///////////////////////////////////////////////////////////////////////////

    const Eigen::SparseVector<double>&
    point_area_gradient(const size_t pi) const
    {
        return m_point_area_jacobian[pi];
    }

    const Eigen::SparseVector<double>& edge_area_gradient(const size_t ei) const
    {
        return m_edge_area_jacobian[ei];
    }
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

    double point_area(const size_t pi) const { return m_point_areas[pi]; }
    const Eigen::VectorXd& point_areas() const { return m_point_areas; }

    double edge_area(const size_t ei) const { return m_edge_areas[ei]; }
    const Eigen::VectorXd& edge_areas() const { return m_edge_areas; }

    ///////////////////////////////////////////////////////////////////////////

    static std::vector<bool> construct_is_on_surface(
        const int num_vertices, const Eigen::MatrixXi& edges);

    /// @brief Helper function that automatically builds include_vertex using construct_is_on_surface.
    /// @param full_vertices_at_rest The full vertices at rest.
    /// @param edges The edge matrix of mesh.
    /// @param faces The face matrix of mesh.
    /// @return Constructed CollisionMesh.
    static CollisionMesh build_from_full_mesh(
        const Eigen::MatrixXd& full_vertices_at_rest,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces)
    {
        return CollisionMesh(
            construct_is_on_surface(full_vertices_at_rest.rows(), edges),
            full_vertices_at_rest, edges, faces);
    }

    /// @brief Construct a matrix that maps from the faces' edges to rows in the edges matrix.
    /// @param faces The face matrix of mesh.
    /// @param edges The edge matrix of mesh.
    /// @return Matrix that maps from the faces' edges to rows in the edges matrix.
    static Eigen::MatrixXi construct_faces_to_edges(
        const Eigen::MatrixXi& faces, const Eigen::MatrixXi& edges);

    /// A function that takes two vertex IDs (row numbers in V) and returns true
    /// if the vertices (and faces or edges containing the vertices) can
    /// collide. By default all primitives can collide with all other
    /// primitives.
    std::function<bool(size_t, size_t)> can_collide = default_can_collide;

protected:
    ///////////////////////////////////////////////////////////////////////////
    // Helper initialization functions

    /// Initialize the selection matrix from full vertices/DOF to collision
    /// vertices/DOF.
    void init_selection_matrices(const int dim);

    /// Initialize point-point and edge-point adjacencies.
    void init_adjacencies();

    /// Initialize point and edge areas.
    void init_areas();

    /// Convert a matrix meant for M_V * V to M_dof * x by duplicating the
    /// entries dim times.
    static Eigen::SparseMatrix<double> vertex_matrix_to_dof_matrix(
        const Eigen::SparseMatrix<double>& M_V, int dim);

    ///////////////////////////////////////////////////////////////////////////

    /// The full vertex positions at rest.
    Eigen::MatrixXd m_full_vertices_at_rest;
    /// The vertex positions at rest.
    Eigen::MatrixXd m_vertices_at_rest;
    /// Edges as rows of indicies into V.
    Eigen::MatrixXi m_edges;
    /// Triangular faces as rows of indicies into V.
    Eigen::MatrixXi m_faces;
    /// Map from F edges to rows of E.
    Eigen::MatrixXi m_faces_to_edges;

    /// Map from full vertices to collision vertices.
    /// @note Negative values indicate full vertex is dropped.
    Eigen::VectorXi m_full_vertex_to_vertex;
    /// Map from collision vertices to full vertices.
    Eigen::VectorXi m_vertex_to_full_vertex;

    /// Selection matrix S ∈ ℝ^{collision×full} for vertices
    Eigen::SparseMatrix<double> m_select_vertices;
    /// Selection matrix S ∈ ℝ^{collision×full} for DOF
    Eigen::SparseMatrix<double> m_select_dof;

    /// Mapping from full displacements to collision displacements
    /// @note this is premultiplied by m_select_vertices
    Eigen::SparseMatrix<double> m_displacement_map;
    /// Mapping from full displacements to collision displacements
    /// @note this is premultiplied by m_select_dof
    Eigen::SparseMatrix<double> m_displacement_dof_map;

    /// Points adjacent to points
    std::vector<unordered_set<int>> m_point_point_adjacencies;
    /// Edges adjacent to edges
    std::vector<unordered_set<int>> m_edge_point_adjacencies;

    // std::vector<std::vector<int>> m_vertices_to_faces;
    // std::vector<std::vector<int>> m_vertices_to_edges;

    /// Is point on the boundary of the triangle mesh in 3D or polyline in 2D?
    std::vector<bool> m_is_point_on_boundary;

    /// Point areas
    /// 2D: 1/2 sum of length of connected edges
    /// 3D: 1/3 sum of area of connected triangles
    Eigen::VectorXd m_point_areas;
    /// Edge areas
    /// 3D: 1/3 sum of area of connected triangles
    Eigen::VectorXd m_edge_areas;

    // Stored as a std::vector so it is easier to access the rows directly.
    /// The rows of the Jacobian of the point areas vector.
    std::vector<Eigen::SparseVector<double>> m_point_area_jacobian;
    /// The rows of the Jacobian of the edge areas vector.
    std::vector<Eigen::SparseVector<double>> m_edge_area_jacobian;

private:
    /// By default all primitives can collide with all other primitives.
    static int default_can_collide(size_t, size_t) { return true; }
};

} // namespace ipc
